"""Compare TF and IR across different mesh resolutions."""

import numpy as np
import scipy.io
import scipy.signal
import scipy.fft
import matplotlib.pyplot as plt
import os

# ---- Config ----
OUTPUT_DIR = r"C:\Masters\DGBABY\edg-acoustics\examples\shoebox\output"
FS_TARGET = 44100
FREQ_LIMIT = 150
ROOM_DIMS = [5, 4, 3]
c0 = 343.0

# Define files and their labels
RUNS = {
    "lc = 0.8": "shoebox_lc08_freq100_2s.mat"
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

    print(f"R={R:.3f}m, halfwidth={halfwidth:.4f}")
    print(f"p_free[0]={p_free[0]:.6f}, prec[0]={prec[0]:.6f}, ratio={prec[0]/p_free[0]:.4f}")
    print(f"p_free peak at t={time_vector[np.argmax(p_free)]*1000:.2f}ms, expected {R/c0*1000:.2f}ms")

    # FFT
    TR_original = scipy.fft.rfft(prec)
    TR_free     = scipy.fft.rfft(p_free)
    freqs       = scipy.fft.rfftfreq(n_samples, dt_sim)

    # Monopole reference
    wavenumber = 2 * np.pi * freqs / c0
    monopole   = np.exp(-1j * wavenumber * R) / (4 * np.pi * R)

    # Regularised deconvolution with frequency mask
    f_max     = 200.0
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

    print(f"IR peak value: {np.max(np.abs(IR_corrected)):.6f}")
    print(f"IR peak time:  {time_vector[np.argmax(np.abs(IR_corrected))]*1000:.2f}ms")

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
    }

# ---- Load all runs ----
results = {}
for label, filename in RUNS.items():
    try:
        results[label] = load_and_process(OUTPUT_DIR, filename, FS_TARGET)
    except FileNotFoundError:
        print(f"WARNING: {filename} not found, skipping.")

for label, res in results.items():
    safe_label = label.replace(" ", "_").replace("=", "").replace(".", "")  # e.g. "lc__08"
    
    # Save as .npz (numpy format)
    np.savez(
        os.path.join(OUTPUT_DIR, f"TR_corrected_lc1_200Hz_2s{safe_label}.npz"),
        TR=res["TR"],
        freqs=res["freqs"],
        IR_resampled=res["IR_resampled"],
        t_resampled=res["t_resampled"]
    )

modes = compute_room_modes(ROOM_DIMS, c0=c0, f_max = FREQ_LIMIT)

# ---- Plot overlaid raw IR ----
plt.figure(figsize=(12, 4))
for label, res in results.items():
    plt.plot(res["t_raw"], res["prec"].squeeze(), label=label, alpha=0.8)
plt.xlabel("Time (s)")
plt.ylabel("Pressure")
plt.title("Raw IR comparison")
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig(os.path.join(OUTPUT_DIR, "compare_raw_IR_lc1_lc1.5_test.png"), dpi=150, bbox_inches="tight")
plt.show()

#plot corrected IR

# ---- Plot overlaid resampled IR ----
plt.figure(figsize=(12, 4))
for label, res in results.items():
    plt.plot(res["t_resampled"], res["IR_resampled"], label=label, alpha=0.8)
plt.xlabel("Time (s)")
plt.ylabel("Pressure")
plt.title(f"Resampled IR comparison @ {FS_TARGET} Hz")
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig(os.path.join(OUTPUT_DIR, "compare_resampled_IR_lc1_lc1.5_test.png"), dpi=150, bbox_inches="tight")
plt.show()


# ---- Plot overlaid Transfer Functions ----
fig, ax = plt.subplots(figsize=(14, 5))
for label, res in results.items():
    ax.plot(
        res["freqs"][res["pos_idx"]],
        20 * np.log10(np.abs(res["TR"][res["pos_idx"]]) + 1e-12),
        label=label, alpha=0.8
    )

for f_mode, nx, ny, nz in modes:
    ax.axvline(x=f_mode, color='red', alpha=0.3, linewidth=0.8, linestyle='--')
    ax.text(f_mode, -58, f"({nx},{ny},{nz})", fontsize=5, rotation=90,
            color='red', ha='center', va='bottom')

ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Magnitude (dB)")
ax.set_title(f"DG Transfer Function with Theoretical Room Modes ({ROOM_DIMS[0]}x{ROOM_DIMS[1]}x{ROOM_DIMS[2]}m)")
ax.set_xlim(0, FREQ_LIMIT)
ax.set_ylim(-60, 10)
ax.legend()
ax.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.show()

print("Done!")

# ---- Plot p_free to verify Gaussian reconstruction ----
plt.figure(figsize=(12, 4))
for label, res in results.items():
    t_raw = res["t_raw"]
    p_free = res["corrected"]["p_free"]
    plt.plot(t_raw, p_free, label=label, alpha=0.8)
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Reconstructed free-field Gaussian pulse (p_free)")
plt.xlim([0, 0.05])  # zoom to first 50ms — pulse should peak at t=R/c0 ≈ 4.3ms
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ---- Plot corrected IR (before resampling) ----
plt.figure(figsize=(12, 4))
for label, res in results.items():
    t_raw = res["t_raw"]
    IR_corrected = res["corrected"]["IR_corrected"]
    plt.plot(t_raw, IR_corrected, label=label, alpha=0.8)
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Corrected IR (before resampling)")
plt.xlim([0, 2.0])  # zoom to first 100ms
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ---- Plot corrected TF (before resampling) ----
plt.figure(figsize=(12, 4))
for label, res in results.items():
    TR_corrected = res["corrected"]["TR_corrected"]
    freqs        = res["corrected"]["freqs"]  # already rfftfreq — use directly
    plt.plot(
        freqs,
        20 * np.log10(np.abs(TR_corrected) + 1e-12),
        label=label, alpha=0.8
    )
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.title("Corrected TF (before resampling, at simulation fs)")
plt.xlim([0, 250])  # only plot up to f_max
plt.ylim(-60, 10)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()



""" to plot the result from main_HF.py
# Load the result
script_dir = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join(script_dir, "result.mat")

print(f"Loading from: {result_path}")
data = scipy.io.loadmat(result_path)

prec = data["prec"]  # shape: [n_rec, n_timesteps]
print(f"prec shape: {prec.shape}")

# You need dt to build the time axis
# dt = CFL * dtscale, but since you didn't save it, 
# either save it too or reconstruct approximately
# For now, load it from the simulation parameters
dt = 0.2 / prec.shape[1]  # total_time / n_timesteps as approximation
t = numpy.arange(prec.shape[1]) * dt

# Plot each receiver
for i in range(prec.shape[0]):
    plt.figure()
    plt.plot(t, prec[i, :])
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure")
    plt.title(f"Raw Pressure at Receiver {i+1}")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    """

"""
#to plot the result from main.py

script_dir = os.path.dirname(os.path.abspath(__file__))
result_path = os.path.join("C:\Masters\DGBABY\edg-acoustics\examples\scenario1\output", "result_meshcoarser.mat")

print(f"Loading from: {result_path}")
data = scipy.io.loadmat(result_path)

IR = data["IR"] #impulse response data
TR = data["TR"] #transfer function
freqs = data["freqs"].squeeze() #freq vector
dt_new = float(data["dt_new"].item()) #upsampled dt
sampling_freq = float(data["sampling_freq"].item()) #should be 44100??
IR_uncorrected = data["IR_Uncorrected"] #raw prec before correction
dt_simulation = float(data["dt_simulation"].item()) #original simulation dt


print(f"Sampling frequency: {sampling_freq} Hz")
print(f"Simulation dt: {dt_simulation:.6f} s  =>  {1/dt_simulation:.1f} Hz")
print(f"IR shape: {IR.shape}")
print(f"Number of receivers: {IR.shape[0]}")
print(f"uncorrected IR shape: {IR_uncorrected.shape}")
print(f"corrected IR shape: {IR.shape}")

#dis shit dont work
#plot corrected IR
t_new = np.arange(IR.shape[1]) * dt_new
plt.plot(t_new, IR[0])
plt.xlabel("time (s)")
plt.ylabel("amplitude")
plt.title("corrected RIR (after source spectrum correction)")
plt.show()


#plot uncorrected (raw) IR (FROM MAIN.PY)
t_old = np.arange(IR_uncorrected.shape[1]) * dt_simulation
plt.plot(t_old, IR_uncorrected[0])
plt.xlabel("time (s)")
plt.ylabel("uncorrected amplitude")
plt.title("uncorrected RIR")
plt.show()

#resample uncorrected IR (FROM MAIN.PY)
import scipy.signal

fs_dg = 1/dt_simulation
fs_target = 44100 

n_samples_new = int(IR_uncorrected[0].shape[0] * fs_target / fs_dg)
IR_dg_resampled = scipy.signal.resample(IR_uncorrected[0], n_samples_new)

print(f"DG IR resampled from {fs_dg: .1f} Hz to {fs_target} Hz")
print(f"Original samples: {IR_uncorrected.shape[1]}, resampled: {n_samples_new}")

dt_resampled = 1 / fs_target
t_resampled = np.arange(len(IR_dg_resampled)) * dt_resampled 

plt.plot(t_resampled, IR_dg_resampled)
plt.title("Resampled uncorrected IR")
plt.xlabel("time (s)")
plt.ylabel("pressure")
plt.show()

#taking FFT of resampled uncorrected IR

n_fft = len(IR_dg_resampled)
TR_uncorrected = scipy.fft.fft(IR_dg_resampled, n=n_fft)
freqs_resampled = scipy.fft.fftfreq(n_fft, dt_resampled)

#plotting the uncorrected, resampled transfer function
pos_idx = freqs_resampled >= 0

plt.figure(figsize=(12, 4))
plt.plot(freqs_resampled[pos_idx], 20 * np.log10(np.abs(TR_uncorrected[pos_idx]) + 1e-12))
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.title("Transfer Function (resampled to 44100 Hz), lc = 1.5")
plt.xlim([0, 300])
plt.grid(True)
plt.tight_layout()
save_path = os.path.join("C:\Masters\DGBABY\edg-acoustics\examples\scenario1\output", "Scenario1_TF_maxfreq200_lc15.png")
plt.savefig(save_path, dpi=150, bbox_inches="tight")
plt.show()
"""


