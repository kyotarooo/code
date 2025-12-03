import meshio
import sys
import os
import glob


def convert_vtk_to_vtu(input_vtk_filename, output_vtu_filename):
    # Read the ASCII VTK file
    mesh = meshio.read(input_vtk_filename)

    # Write the VTU file in binary format
    meshio.write(output_vtu_filename, mesh, file_format="vtu", binary=True)

    print(f"Converted {input_vtk_filename} → {output_vtu_filename} (binary VTU)")


if __name__ == "__main__":
    vtk_files = glob.glob("elem_*.vtk")
    for vtk_file in vtk_files:
        if os.path.isfile(vtk_file):
            input_file = vtk_file
            output_file = vtk_file.removesuffix(".vtk") + ".vtu"
            convert_vtk_to_vtu(input_file, output_file)
            os.remove(vtk_file)
    
    file_name = "time.out"
    if os.path.isfile(file_name):
        with open(file_name, "r") as f:
            lines = f.readlines()
        file_name = "elem.pvd"
        with open(file_name, "w") as f:
            f.write("<VTKFile type=\"Collection\">\n")
            f.write("  <Collection>\n")
            for line in lines:
                id = int(line.split()[0])
                time = float(line.split()[1])
                f.write(f"    <DataSet part=\"0\" timestep=\"{time}\" file=\"elem_{id:04d}.vtu\" />\n")
            f.write("  </Collection>\n")
            f.write("</VTKFile>\n")
