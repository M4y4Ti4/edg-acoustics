#main code for room with partial wall 

""" This is the main script for the shoebox"""

# region Import Libraries
import os
import glob
import numpy
import scipy.io
import edg_acoustics
import time


# endregion

# --------------------
# Block 1: User input
# --------------------
rho0 = 1.213  # density of air at 20 degrees Celsius in kg/m^3
c0 = 343  # speed of sound in air at 20 degrees Celsius in m/s

BC_labels = {
    "carpet": 2,
    "ceiling": 3,
}  # predefined labels for boundary conditions. please assign an arbitrary int number to each type of boundary condition, e.g. hard wall, carpet, panel. The number should be unique for each type of boundary condition and should match the physical surface number in the .geo mesh file. The string should be the same as the material name in the .mat file (at least for the first few letters).

real_valued_impedance_boundary = [
    {"label":  2, "RI": 0.9849},
    {"label": 3, "RI": 0.9592},
]# extra labels for real-valued impedance boundary condition, if needed. The label should be the similar to the label in BC_labels. Since it's frequency-independent, only "RI", the real-valued reflection coefficient, is required. If not needed, just clear the elements of this list and keep the empty list.

mesh_used = "room_with_wall_fixed"
mesh_name = f"{mesh_used}.msh"  # name of the mesh file. The mesh file should be in the same folder as this script.
monopole_xyz = numpy.array([-1.04, 2.5, 1.62])  # x,y,z coordinate of the source in the room
freq_upper_limit = 100  # upper limit of the frequency content of the source signal in Hz. The source signal is a Gaussian pulse with a frequency content up to this limit.

# Approximation degrees
Nx = 4  # in space
Nt = 4  # in time
CFL = 0.5  # CFL number, default is 0.5.
recx = numpy.array([-4.5])
recy = numpy.array([2.5])
recz = numpy.array([1.62])
rec = numpy.vstack((recx, recy, recz))  # dim:[3,n_rec]

impulse_length = 2.0  # total simulation time in seconds
save_every_Nstep = 1  # save the results every N steps
temporary_save_Nstep = 500  # save the results every N steps temporarily during the simulation. The temporary results will be saved in the root directory of this repo.

#define output directory
output_dir = os.path.join(os.path.split(os.path.abspath(__file__))[0], "output")
os.makedirs(output_dir, exist_ok=True)  # creates folder if it doesn't exist
result_filename = "shoebox_lc25_freq200_2s"  # name of the result file. The result file will be saved in the same folder as this script. The result file will be saved in .mat format.

# --------------------------------------------------------------------------------
# Block 2: Initialize the simulation，run the simulation and save the results
# --------------------------------------------------------------------------------

# load Boundary conditions and parameters
BC_para = real_valued_impedance_boundary

# mesh_data_folder is the current folder by default
mesh_data_folder = os.path.split(os.path.abspath(__file__))[0]
mesh_filename = os.path.join(mesh_data_folder, mesh_name)
mesh = edg_acoustics.Mesh(mesh_filename, BC_labels)


IC = edg_acoustics.Monopole_IC(monopole_xyz, freq_upper_limit)

sim = edg_acoustics.AcousticsSimulation(rho0, c0, Nx, mesh, BC_labels)

flux = edg_acoustics.UpwindFlux(rho0, c0, sim.n_xyz)
AbBC = edg_acoustics.AbsorbBC(sim.BCnode, BC_para)

sim.init_BC(AbBC)
sim.init_IC(IC)
sim.init_Flux(flux)
sim.init_rec(
    rec, "scipy"
)  # brute_force or scipy(default) approach to locate the receiver points in the mesh

simulation_start = time.time() #tracking the time of the simulation

tsi_time_integrator = edg_acoustics.TSI_TI(sim.RHS_operator, sim.dtscale, CFL, Nt=3)
sim.init_TimeIntegrator(tsi_time_integrator)
sim.time_integration(
    total_time=impulse_length,
    delta_step=save_every_Nstep,
    save_step=temporary_save_Nstep,
    format="mat",
)

simulation_elapsed = time.time() - simulation_start

hours = int(simulation_elapsed // 3600)
minutes = int((simulation_elapsed % 3600) // 60)
seconds = int(simulation_elapsed % 60)
print(f"Time integration: {hours}h {minutes}m {seconds}s")

#saving the raw results directly - no postprocessing
result_path = os.path.join(output_dir, result_filename)
scipy.io.savemat(f"{result_path}.mat", {
    "prec": sim.prec,
    "dt": sim.time_integrator.dt,
    "runtime_seconds": simulation_elapsed,
    "runtime_string": f"{seconds}s",
    "source_xyz": sim.IC.source_xyz,
    "halfwidth": sim.IC.halfwidth,
    "rec": rec,
    "c0": c0,
    "rho0": rho0,
    "mesh_name": mesh_name,
    "Nx": Nx,
    "CFL": CFL,
    "N_tets": sim.N_tets,
    "impulse_length": impulse_length,
})

print(f"Results saved to: {result_path}.mat")  # this was also wrong before
print("Finished!")
