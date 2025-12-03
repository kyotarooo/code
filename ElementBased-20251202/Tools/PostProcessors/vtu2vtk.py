import meshio
import sys
import os
import glob


def convert_vtu_to_vtk(input_vtu_filename, output_vtk_filename):
    # Read the VTU file (binary or ASCII)
    mesh = meshio.read(input_vtu_filename)

    # Write the legacy VTK file in ASCII format
    meshio.write(
        output_vtk_filename,
        mesh,
        file_format="vtk",  # legacy VTK format
        binary=False,  # ASCII output
    )

    print(f"Converted {input_vtu_filename} → {output_vtk_filename} (ASCII VTK)")


if __name__ == "__main__":
    vtu_files = glob.glob("elem_*.vtu")
    for vtu_file in vtu_files:
        if os.path.isfile(vtu_file):
            input_file = vtu_file
            output_file = vtu_file.removesuffix(".vtu") + ".vtk"
            print(input_file)
            print(output_file)
            convert_vtu_to_vtk(input_file, output_file)
            os.remove(vtu_file)
