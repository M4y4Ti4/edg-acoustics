#convert geo file to msh file 

"""
import gmsh
gmsh.initialize()
gmsh.open(r"C:\Masters\Hybrid\DGsim\examples\wall\room_with_wall.geo")
gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(r"C:\Masters\Hybrid\DGsim\examples\wall\room_with_wall.msh")
gmsh.finalize()
print("Done!")
"""
"""
import gmsh
import meshio

gmsh.initialize()
gmsh.model.add("room_with_wall")

# Use OpenCASCADE kernel
factory = gmsh.model.occ

# Room outer box: corner at (-6.68, 0, 0), dimensions (6.68, 5.98, 4.06)
factory.addBox(-6.68, 0.0, 0.0, 6.68, 5.98, 4.06)  # tag=1

# Partial wall box: thin wall at x=-3.926 to -3.706, y=0 to 2.72, full height
factory.addBox(-3.926297, 0.0, 0.0, 0.22, 2.72, 4.06)  # tag=2

# Synchronise before Boolean operation
factory.synchronize()

# Cut wall from room — removes wall volume, keeps room with cavity
out, _ = factory.cut([(3, 1)], [(3, 2)], removeTool=True)

factory.synchronize()

# Check result
volumes  = gmsh.model.getEntities(3)
surfaces = gmsh.model.getEntities(2)
print(f"Volumes after cut:  {volumes}")
print(f"Surfaces after cut: {surfaces}")

# Auto-identify surfaces by position
floor_tags   = []
ceiling_tags = []

for dim, tag in surfaces:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    zc = (zmin + zmax) / 2
    yc = (ymin + ymax) / 2
    
    if zmax < 0.01:                    # floor: z ≈ 0
        floor_tags.append(tag)
    else:
        ceiling_tags.append(tag)       # everything else

print(f"Floor surfaces:   {floor_tags}")
print(f"Ceiling surfaces: {ceiling_tags}")

# Set mesh size
lc = 0.5
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)

# Generate 3D mesh
gmsh.option.setNumber("Mesh.Algorithm",   1)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)
gmsh.option.setNumber("Mesh.Optimize",    1)

gmsh.model.mesh.generate(3)

# Assign named physical groups matching BC_labels
vol_tags = [v[1] for v in volumes]
gmsh.model.addPhysicalGroup(3, vol_tags,     name="air")
gmsh.model.addPhysicalGroup(2, floor_tags,   name="carpet")
gmsh.model.addPhysicalGroup(2, ceiling_tags, name="ceiling")

gmsh.write("room_with_wall.msh")
gmsh.finalize()

# Verify
m = meshio.read("room_with_wall.msh")
print("\nCell types:", list(m.cells_dict.keys()))
for k, v in m.cells_dict.items():
    print(f"  {k}: {len(v)} elements")
"""
"""
Fix MSH file dtype issue that causes scipy Delaunay to fail.
Run this after generating room_with_wall.msh with the OpenCASCADE script.
It rewrites the mesh with float64 node coordinates and int32 element indices.
"""
import meshio
import numpy as np
 
input_path  = r"room_with_wall.msh"
output_path = r"room_with_wall_fixed.msh"
 
print(f"Reading {input_path}...")
mesh = meshio.read(input_path)
 
# Force all point coordinates to float64
mesh.points = mesh.points.astype(np.float64)
 
# Force all cell connectivity arrays to int32
new_cells = []
for cell_block in mesh.cells:
    new_data = cell_block.data.astype(np.int32)
    new_cells.append(meshio.CellBlock(cell_block.type, new_data))
mesh.cells = new_cells
 
# Also fix cell_data if present
new_cell_data = {}
for key, val_list in mesh.cell_data.items():
    new_cell_data[key] = [v.astype(np.int32) for v in val_list]
mesh.cell_data = new_cell_data
 
print(f"Points dtype:  {mesh.points.dtype}")
print(f"Cell dtypes:   {[c.data.dtype for c in mesh.cells]}")
print(f"Cell types:    {[c.type for c in mesh.cells]}")
print(f"Num tets:      {next(c.data.shape[0] for c in mesh.cells if c.type == 'tetra')}")
 
# Write fixed mesh in MSH 2.2 format
meshio.write(output_path, mesh, file_format="gmsh22")
print(f"Fixed mesh written to {output_path}")
print("Update mesh_name in wall_main.py to 'room_with_wall_fixed.msh'")