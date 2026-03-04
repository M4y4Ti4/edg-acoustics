import numpy
import scipy.io
import matplotlib.pyplot as plt
import os

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