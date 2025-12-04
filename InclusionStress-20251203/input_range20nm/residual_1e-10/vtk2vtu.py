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
    vtk_files = glob.glob("section.vtk")
    for vtk_file in vtk_files:
        if os.path.isfile(vtk_file):
            input_file = vtk_file
            output_file = vtk_file.removesuffix(".vtk") + ".vtu"
            convert_vtk_to_vtu(input_file, output_file)
            os.remove(vtk_file)
