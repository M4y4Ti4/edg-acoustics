import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import os
import scipy.fft

#plotting from raw results (without postprocessing)
output_dir = r"C:\Masters\DGBABY\edg-acoustics\examples\scenario1\output"
result_filename = "scenario1_coarser.mat"

#loading data
result_path = result_path = os.path.join(output_dir, result_filename)
data = scipy.io.loadmat(result_path)

prec = data["prec"]
dt_sim = data["dt"].item()
fs_dg = 1 / dt_sim
fs_target = 44100

print(f"Simulation dt: {dt_sim:.6f} s  =>  {fs_dg:.1f} Hz")
print(f"prec shape: {prec.shape}")
print(f"Duration: {prec.shape[1] * dt_sim:.2f} s")
print(f"Runtime: {data['runtime_string']}")

#plotting the unsampled, uncorrected RIR
t_old = np.arange(prec.shape[1]) * dt_sim
plt.plot(t_old, prec[0])
plt.xlabel("Time (s)")
plt.ylabel("Pressure")
plt.title(f"Uncorrected, not resampled IR, lc = 1.5, max freq = 200")
plt.show()

#resample IR
n_samples_new = int(prec.shape[1] * fs_target / fs_dg)
IR_resampled = scipy.signal.resample(prec[0], n_samples_new)
dt_resampled = 1 / fs_target
t_resampled = np.arange(len(IR_resampled)) * dt_resampled

#plot resampled, uncorrected IR
plt.plot(t_resampled, IR_resampled)
plt.xlabel("Time (s)")
plt.ylabel("Pressure")
plt.title(f"Resampled IR to 44100 Hz, lc = 1.5")
plt.grid()
plt.show()

#apply fourier transform
n_fft = len(IR_resampled)
TR = scipy.fft.fft(IR_resampled, n=n_fft)
freqs = scipy.fft.fftfreq(n_fft, dt_resampled)
pos_idx = freqs >= 0

#plotting transfer function
plt.plot(freqs[pos_idx], 20 * np.log10(np.abs(TR[pos_idx]) + 1e-12))
plt.xlabel("freq")
plt.ylabel("mag")
plt.xlim([0, 300])
plt.title("TF")
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


