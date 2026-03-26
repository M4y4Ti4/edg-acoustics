"""Compare TF and IR across different mesh resolutions."""

import numpy as np
import scipy.io
import scipy.signal
import scipy.fft
import matplotlib.pyplot as plt
import os

# ---- Config ----
OUTPUT_DIR = r"C:\Masters\DGsim\edg-acoustics\examples\shoebox\output"
FS_TARGET = 44100
FREQ_LIMIT = 200
ROOM_DIMS = [5, 4, 3]
c0 = 343.0

# Define files and their labels
RUNS = {
    "lc = 1.0": "shoebox_lc1_freq200_2s_highabs.mat"
}

def compute_room_modes(room_dims, c0=343.0, f_max=150, n_modes = 5):
    lx, ly, lz = room_dims
    modes = []
    for nx in range(n_modes):
        for ny in range(n_modes):
            for nz in range(n_modes):
                if nx == ny == nz == 0: 
                    continue
                if sum([nx > 0, ny > 0, nz > 0]) > 1:
                    continue
                f_mode = (c0 / 2) * np.sqrt((nx/lx)**2 + (ny/ly)**2 + (nz/lz)**2)
                if f_mode <= f_max: 
                    modes.append((f_mode, nx, ny, nz))
    modes.sort()
    return modes

def apply_correction(prec, dt_sim, source_xyz, rec_xyz, c0=343.0, halfwidth=None):

    prec       = np.array(prec).squeeze()
    source_xyz = np.array(source_xyz).squeeze()
    rec_xyz    = np.array(rec_xyz).squeeze()
    halfwidth  = float(np.array(halfwidth).squeeze())
    

    n_samples   = len(prec)
    time_vector = np.arange(n_samples) * dt_sim
    fs          = 1.0 / dt_sim
    R           = np.sqrt(np.sum((source_xyz - rec_xyz)**2))

    # Free-field solution for pure Gaussian pressure IC, zero velocity IC
    p_free = 0.5 * (
        ((R + c0 * time_vector) / R) *
        np.exp(-np.log(2) * (R - c0 * time_vector)**2 / halfwidth**2) +
        ((R - c0 * time_vector) / R) *
        np.exp(-np.log(2) * (R + c0 * time_vector)**2 / halfwidth**2)
    )

    # FFT
    TR_original = scipy.fft.rfft(prec)
    TR_free     = scipy.fft.rfft(p_free)
    freqs       = scipy.fft.rfftfreq(n_samples, dt_sim)

    # Monopole reference
    wavenumber = 2 * np.pi * freqs / c0
    monopole   = np.exp(-1j * wavenumber * R) / (4 * np.pi * R)

    # Regularised deconvolution with frequency mask
    f_max     = 300.0
    epsilon   = 1e-4 * np.max(np.abs(TR_free)**2)
    freq_mask = freqs <= f_max

    TR_corrected = np.zeros(len(TR_original), dtype=complex)
    TR_corrected[freq_mask] = (
        TR_original[freq_mask] * np.conj(TR_free[freq_mask]) /
        (np.abs(TR_free[freq_mask])**2 + epsilon)
    ) * monopole[freq_mask]

    IR_corrected = scipy.fft.irfft(TR_corrected, n=n_samples)

    # Window to suppress noise tail
    window = np.ones(n_samples)
    t_start = int(0.70 * fs)
    t_end   = int(0.80 * fs)
    if t_end < n_samples:
        taper = np.arange(t_end - t_start)
        window[t_start:t_end] = 0.5 * (1 + np.cos(np.pi * taper / (t_end - t_start)))
        window[t_end:] = 0.0
    IR_corrected *= window

    

    return {
        "IR_corrected": IR_corrected,
        "TR_corrected": TR_corrected,
        "TR_original":  TR_original,
        "TR_free":      TR_free,
        "freqs":        freqs,
        "p_free":       p_free,
    }

def load_and_process(output_dir, filename, fs_target):
    """Load, resample and compute TF for a single result file."""
    path = os.path.join(output_dir, filename)
    data = scipy.io.loadmat(path)
    prec = data["prec"].squeeze()
    dt_sim = data["dt"].item()
    t_raw = np.arange(len(prec)) * dt_sim

    peak_idx = np.argmax(prec)
    print(f"Raw prec peak value: {prec[peak_idx]:.6f}")
    print(f"Raw prec peak time:  {t_raw[peak_idx]*1000:.2f}ms")
    print(f"prec[0]:             {prec[0]:.6f}")
    print(f"Expected direct arrival: {1.476/343*1000:.2f}ms")
    print("All mat keys and values:")
    for k, v in data.items():
        if not k.startswith('__'):
            print(f"  {k}: {np.array(v).squeeze()}")

    prec = data["prec"].squeeze()
    dt_sim = data["dt"].item()
    fs_dg = 1 / dt_sim
    source_pos = data["source_xyz"].squeeze()
    rec_pos = data["rec"].squeeze()
    halfwidth = float(data["halfwidth"].squeeze())
    print(halfwidth)

    corrected = apply_correction(prec, dt_sim, source_xyz = source_pos, rec_xyz = rec_pos, halfwidth = halfwidth)
        
    IR_corrected = corrected["IR_corrected"]
    TR_corrected = corrected["TR_corrected"]

    # Time axis (raw)
    t_raw = np.arange(len(prec)) * dt_sim

    # Resample
    n_samples_new = int(len(IR_corrected) * fs_target / fs_dg)
    IR_resampled = scipy.signal.resample(IR_corrected, n_samples_new)
    dt_resampled = 1 / fs_target
    t_resampled = np.arange(len(IR_resampled)) * dt_resampled

    # Transfer function
    n_fft = len(IR_resampled)
    TR = scipy.fft.fft(IR_resampled, n=n_fft)
    freqs = scipy.fft.fftfreq(n_fft, dt_resampled)
    pos_idx = freqs >= 0

    return {
    "t_raw": t_raw,
    "prec": prec,
    "t_resampled": t_resampled,
    "IR_resampled": IR_resampled,
    "TR": TR,
    "freqs": freqs,
    "pos_idx": pos_idx,
    "corrected": corrected,
    "TR_corrected": TR_corrected
}


def post_process_output(prec, dt_sim, source_xyz, rec_xyz, halfwidth, fs_target=44100):
    """Same as load_and_process but accepts data directly instead of loading from file"""
    import numpy 
    fs_dg = 1 / dt_sim
    
    corrected = apply_correction(prec, dt_sim, source_xyz=source_xyz, rec_xyz=rec_xyz, halfwidth=halfwidth)
    
    IR_corrected = corrected["IR_corrected"]
    
    # Resample
    n_samples_new = int(len(IR_corrected) * fs_target / fs_dg)
    IR_resampled  = scipy.signal.resample(IR_corrected, n_samples_new)
    t_resampled   = numpy.arange(len(IR_resampled)) / fs_target

    n_fft = len(IR_resampled)
    TF = scipy.fft.fft(IR_resampled, n=n_fft)
    freqs = scipy.fft.rfftfreq(n_fft, 1/fs_target)

    return {
        "IR": IR_resampled,
        "TF": TF,
        "t_resampled": t_resampled,
        "freqs": freqs,
        "fs": fs_target
    }


