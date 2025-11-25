# ======================================================
# inputgen.py: Generator of input files for incstrgen
# ------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# ======================================================
import TkEasyGUI as eg


def write_grid(values, selected_folder):
    with open(f"{selected_folder}/grid.inp", "w") as f:
        n = [int(values["-NX-"]), int(values["-NY-"]), int(values["-NZ-"])]
        f.write(" ".join(map(str, n)) + "\n")
        n = [int(values["-SUBNX-"]), int(values["-SUBNY-"]), int(values["-SUBNZ-"])]
        f.write(" ".join(map(str, n)) + "\n")


def write_material(values, selected_folder):
    with open(f"{selected_folder}/material.inp", "w") as f:
        tensor = [
            [float(values["-OXX-"]), float(values["-OXY-"]), float(values["-OXZ-"])],
            [float(values["-OYX-"]), float(values["-OYY-"]), float(values["-OYZ-"])],
            [float(values["-OZX-"]), float(values["-OZY-"]), float(values["-OZZ-"])],
        ]
        for vector in tensor:
            f.write(" ".join(map(str, vector)) + "\n")

        periodic = [int(values["-PBCX-"]), int(values["-PBCY-"]), int(values["-PBCZ-"])]
        f.write(" ".join(map(str, periodic)) + "\n")

        size = [float(values["-SX-"]), float(values["-SY-"]), float(values["-SZ-"])]
        f.write(" ".join(map(str, size)) + "\n")

        f.write(f"{values["-SHEAR-"]} {values["-POISSON-"]}\n")

def open_window():
    layout = [
        [eg.Frame("Grid Information", [
        [
            eg.Text("Number of grids"),
            eg.Text("X"),
            eg.Input("10", key="-NX-"),
            eg.Text("Y"),
            eg.Input("10", key="-NY-"),
            eg.Text("Z"),
            eg.Input("10", key="-NZ-"),
        ],
        [
            eg.Text("Number of subgrids"),
            eg.Text("X"),
            eg.Input("0", key="-SUBNX-"),
            eg.Text("Y"),
            eg.Input("0", key="-SUBNY-"),
            eg.Text("Z"),
            eg.Input("0", key="-SUBNZ-"),
        ]],expand_x=True)],
        [eg.Frame("Material Information", [
        [eg.Text("Crystal orientation of volume")],
        [
            eg.Text("XX"),
            eg.Input("1.0", key="-OXX-"),
            eg.Text("XY"),
            eg.Input("0.0", key="-OXY-"),
            eg.Text("XZ"),
            eg.Input("0.0", key="-OXZ-"),
        ],
        [
            eg.Text("YX"),
            eg.Input("0.0", key="-OYX-"),
            eg.Text("YY"),
            eg.Input("1.0", key="-OYY-"),
            eg.Text("YZ"),
            eg.Input("0.0", key="-OYZ-"),
        ],
        [
            eg.Text("ZX"),
            eg.Input("0.0", key="-OZX-"),
            eg.Text("ZY"),
            eg.Input("0.0", key="-OZY-"),
            eg.Text("ZZ"),
            eg.Input("1.0", key="-OZZ-"),
        ],
        [
            eg.Text("Periodic boundary condition"),
            eg.Checkbox("X", True, key="-PBCX-"),
            eg.Checkbox("Y", True, key="-PBCY-"),
            eg.Checkbox("Z", True, key="-PBCZ-"),
        ],
        [
            eg.Text("Simulation volume size"),
            eg.Text("X"),
            eg.Input("0.0", key="-SX-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-SY-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-SZ-"),
        ],
        [eg.Text("Elastic constants")],
        [
            eg.Text("Shear modulus (Pa)"),
            eg.Input("80.0e+09", key="-SHEAR-"),
            eg.Text("Poisson's ratio"),
            eg.Input("0.3", key="-POISSON-"),
        ]],expand_x=True)],
        [eg.Button("Save"), eg.Button("Cancel")],
    ]
    return eg.Window("Input File Generator for incstrgen", layout=layout)


def main():
    window = open_window()
    
    while window.is_alive():
        event, values = window.read()

        if event == "Save":
            window.close()
            selected_folder = eg.popup_get_folder("Select folder")
            if selected_folder:
                write_grid(values, selected_folder)
                write_material(values, selected_folder)
                eg.popup("Input files are saved.")
            else:
                eg.popup_warning("Input files not saved.")
            break
        if event == "Cancel":
            window.close()
            eg.popup_warning("Input files not saved.")
            break


if __name__ == "__main__":
    main()
