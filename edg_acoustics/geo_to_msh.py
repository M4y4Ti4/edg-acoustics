#convert geo file to msh file 
import gmsh
gmsh.initialize()
gmsh.open("C:\Masters\DGBABY\edg-acoustics\examples\scenario1\scenario1_lc1.geo")
gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(r"C:\Masters\DGBABY\edg-acoustics\examples\scenario1\scenario1_lc1.msh")
gmsh.finalize()
print("Done!")