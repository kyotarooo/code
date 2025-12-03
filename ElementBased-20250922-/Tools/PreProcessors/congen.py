# ===================================================================
# congen-gui.py: an interactive condition.inp file generator with GUI
# -------------------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# ===================================================================
import TkEasyGUI as eg

# Default parameters
defaults = {
    "NUM_STEPS": "10000",
    "UMAX": "1.0e-09",
    "DTMIN": "1.0e-12",
    "DTMAX": "1.0e-06",
    "OUTPUT_INTERVAL": "100",
    "ELASTIC_SHEAR_MODULUS": "80.0e+09",
    "POISSONS_RATIO": "0.3",
    "MOBILITY": "1000.0",
    "CORE_RADIUS": "12.5e-10",
    "VOLUME_X": "1.0e-06",
    "VOLUME_Y": "1.0e-06",
    "VOLUME_Z": "1.0e-06",
    "PBC_X": 0,
    "PBC_Y": 0,
    "PBC_Z": 0,
    "REACTION_RADIUS": "15.0e-10",
    "STRESS_XX": "0.0e+06",
    "STRESS_YY": "0.0e+06",
    "STRESS_ZZ": "0.0e+06",
    "STRESS_XY": "0.0e+06",
    "STRESS_YZ": "0.0e+06",
    "STRESS_XZ": "0.0e+06",
    "INDEX_I": "-1",
    "INDEX_J": "-1",
    "STRAIN_RATE": "10000.0",
    "MIN_ELEMENT_LENGTH": "1.0e-08",
    "MAX_ELEMENT_LENGTH": "1.0e-07"
}


# Load existing file or use defaults
def load_existing_file(file_name):
    try:
        with open(file_name, "r") as f:
            data_lines = f.readlines()
        values = {
            "NUM_STEPS": data_lines[0].split()[0],
            "UMAX": data_lines[1].split()[0],
            "DTMIN": data_lines[2].split()[0],
            "DTMAX": data_lines[2].split()[1],
            "OUTPUT_INTERVAL": data_lines[3].split()[0],
            "ELASTIC_SHEAR_MODULUS": data_lines[4].split()[0],
            "POISSONS_RATIO": data_lines[4].split()[1],
            "MOBILITY": data_lines[5].split()[0],
            "CORE_RADIUS": data_lines[6].split()[0],
            "VOLUME_X": data_lines[7].split()[0],
            "VOLUME_Y": data_lines[7].split()[1],
            "VOLUME_Z": data_lines[7].split()[2],
            "PBC_X": int(data_lines[8].split()[0]),
            "PBC_Y": int(data_lines[8].split()[1]),
            "PBC_Z": int(data_lines[8].split()[2]),
            "REACTION_RADIUS": data_lines[9].split()[0],
            "STRESS_XX": data_lines[10].split()[0],
            "STRESS_YY": data_lines[11].split()[1],
            "STRESS_ZZ": data_lines[12].split()[2],
            "STRESS_XY": data_lines[10].split()[1],
            "STRESS_YZ": data_lines[11].split()[2],
            "STRESS_XZ": data_lines[10].split()[2],
            "INDEX_I": data_lines[13].split()[0],
            "INDEX_J": data_lines[13].split()[1],
            "STRAIN_RATE": data_lines[13].split()[2],
            "MIN_ELEMENT_LENGTH": data_lines[14].split()[0],
            "MAX_ELEMENT_LENGTH": data_lines[14].split()[1],
        }
        return values
    except FileNotFoundError:
        eg.popup_warning("No existing condition.inp file found. Using default values.")
        return defaults

def setup_layout(values):
    layout = [
        [
            eg.Frame(
                "Time Step Configurations",
                [
                    [
                        eg.Text("Number of time steps"),
                        eg.Input(values["NUM_STEPS"], key="-NUM_STEPS-"),
                    ],
                    [
                        eg.Text("Maximum displacement for each time step (m)"),
                        eg.Input(values["UMAX"], key="-UMAX-"),
                    ],
                    [
                        eg.Text("Minimum time step size (s)"),
                        eg.Input(values["DTMIN"], key="-DTMIN-"),
                    ],
                    [
                        eg.Text("Maximum time step size (s)"),
                        eg.Input(values["DTMAX"], key="-DTMAX-"),
                    ],
                    [
                        eg.Text("Output interval"),
                        eg.Input(values["OUTPUT_INTERVAL"], key="-OUTPUT_INTERVAL-"),
                    ],
                ],
                expand_x=True,
            )
        ],
        [
            eg.Frame(
                "Material Properties",
                [
                    [
                        eg.Text("Elastic shear modulus (Pa)"),
                        eg.Input(values["ELASTIC_SHEAR_MODULUS"], key="-SHEAR_MODULUS-"),
                    ],
                    [
                        eg.Text("Poisson's ratio"),
                        eg.Input(values["POISSONS_RATIO"], key="-POISSONS_RATIO-"),
                    ],
                    [
                        eg.Text("Mobility (/Ps s)"),
                        eg.Input(values["MOBILITY"], key="-MOBILITY-"),
                    ],
                    [
                        eg.Text("Dislocation core radius (m)"),
                        eg.Input(values["CORE_RADIUS"], key="-CORE_RADIUS-"),
                    ],
                ],
                expand_x=True,
            )
        ],
        [
            eg.Frame(
                "Boundary Conditions",
                [
                    [eg.Text("Simulation volume size (m)")],
                    [
                        eg.Text("X"),
                        eg.Input(values["VOLUME_X"], key="-VOLUME_X-"),
                        eg.Text("Y"),
                        eg.Input(values["VOLUME_Y"], key="-VOLUME_Y-"),
                        eg.Text("Z"),
                        eg.Input(values["VOLUME_Z"], key="-VOLUME_Z-"),
                    ],
                    [
                        eg.Text("Periodic boundary condition"),
                        eg.Checkbox("X", bool(values["PBC_X"]), key="-PBC_X-"),
                        eg.Checkbox("Y", bool(values["PBC_Y"]), key="-PBC_Y-"),
                        eg.Checkbox("Z", bool(values["PBC_Z"]), key="-PBC_Z-"),
                    ],
                    [
                        eg.Text("Direct reaction radius (m)"),
                        eg.Input(values["REACTION_RADIUS"], key="-REACTION_RADIUS-"),
                    ],
                    [eg.Text("Applied stress (Pa)")],
                    [
                        eg.Text("XX"),
                        eg.Input(values["STRESS_XX"], key="-STRESS_XX-"),
                        eg.Text("YY"),
                        eg.Input(values["STRESS_YY"], key="-STRESS_YY-"),
                        eg.Text("ZZ"),
                        eg.Input(values["STRESS_ZZ"], key="-STRESS_ZZ-"),
                    ],
                    [
                        eg.Text("XY"),
                        eg.Input(values["STRESS_XY"], key="-STRESS_XY-"),
                        eg.Text("YZ"),
                        eg.Input(values["STRESS_YZ"], key="-STRESS_YZ-"),
                        eg.Text("XZ"),
                        eg.Input(values["STRESS_XZ"], key="-STRESS_XZ-"),
                    ],
                    [eg.Text("Strain rate test condition")],
                    [
                        eg.Text("Index (i)"),
                        eg.Input(values["INDEX_I"], key="-INDEX_I-", width=5),
                        eg.Text("Index (j)"),
                        eg.Input(values["INDEX_J"], key="-INDEX_J-", width=5),
                        eg.Text("Strain rate (/s)"),
                        eg.Input(values["STRAIN_RATE"], key="-STRAIN_RATE-"),
                    ],
                ],
                expand_x=True,
            )
        ],
        [
            eg.Frame(
                "Simulation Parameters",
                [
                    [
                        eg.Text("Standard element length (m)"),
                        eg.Text("Minimum"),
                        eg.Input(values["MIN_ELEMENT_LENGTH"], key="-MIN_ELEMENT_LENGTH-"),
                        eg.Text("Maximum"),
                        eg.Input(values["MAX_ELEMENT_LENGTH"], key="-MAX_ELEMENT_LENGTH-"),
                    ]
                ],
                expand_x=True,
            )
        ],
        [eg.Button("Save"), eg.Button("Load"), eg.Button("Cancel")],
    ]
    return layout

layout = setup_layout(defaults)
window = eg.Window(
    "Condition.inp file Generator for Element-based PDD Simulations", layout=layout
)

while window.is_alive():
    event, values = window.read()
    if event == "Save":
        window.close()
        selected_dir = eg.popup_get_folder("Select folder")
        with open(f"{selected_dir}/condition.inp", "w") as f:
            f.write(f"{values["-NUM_STEPS-"]}\n")
            f.write(f"{values["-UMAX-"]}\n")
            f.write(f"{values["-DTMIN-"]} {values["-DTMAX-"]}\n")
            f.write(f"{values["-OUTPUT_INTERVAL-"]}\n")
            f.write(f"{values["-SHEAR_MODULUS-"]} {values["-POISSONS_RATIO-"]}\n")
            f.write(f"{values["-MOBILITY-"]}\n")
            f.write(f"{values["-CORE_RADIUS-"]}\n")
            f.write(
                f"{values["-VOLUME_X-"]} {values["-VOLUME_Y-"]} {values["-VOLUME_Z-"]}\n"
            )
            f.write(
                f"{int(values["-PBC_X-"])} {int(values["-PBC_Y-"])} {int(values["-PBC_Z-"])}\n"
            )
            f.write(f"{values["-REACTION_RADIUS-"]}\n")
            f.write(
                f"{values["-STRESS_XX-"]} {values["-STRESS_XY-"]} {values["-STRESS_XZ-"]}\n"
            )
            f.write(
                f"{values["-STRESS_XY-"]} {values["-STRESS_YY-"]} {values["-STRESS_YZ-"]}\n"
            )
            f.write(
                f"{values["-STRESS_XZ-"]} {values["-STRESS_YZ-"]} {values["-STRESS_ZZ-"]}\n"
            )
            f.write(
                f"{values["-INDEX_I-"]} {values["-INDEX_J-"]} {values["-STRAIN_RATE-"]}\n"
            )
            f.write(f"{values["-MIN_ELEMENT_LENGTH-"]} {values["-MAX_ELEMENT_LENGTH-"]}\n")
        with open(f"{selected_dir}/volume.vtk", "w") as f:
            x0 = 0.0
            x1 = float(values["-VOLUME_X-"])
            y0 = 0.0
            y1 = float(values["-VOLUME_Y-"])
            z0 = 0.0
            z1 = float(values["-VOLUME_Z-"])
            f.write("# vtk DataFile Version 3.0\n")
            f.write("0.0\n")
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")
            f.write(f"POINTS 8 float\n")
            f.write(f"{x0:e} {y0:e} {z0:e}\n")
            f.write(f"{x1:e} {y0:e} {z0:e}\n")
            f.write(f"{x1:e} {y1:e} {z0:e}\n")
            f.write(f"{x0:e} {y1:e} {z0:e}\n")
            f.write(f"{x0:e} {y0:e} {z1:e}\n")
            f.write(f"{x1:e} {y0:e} {z1:e}\n")
            f.write(f"{x1:e} {y1:e} {z1:e}\n")
            f.write(f"{x0:e} {y1:e} {z1:e}\n")
            f.write("CELLS 1 9\n")
            f.write("8 0 1 2 3 4 5 6 7\n")
        eg.popup("condition.inp file is generated.")
        break
    if event == "Load":
        window.close()
        selected_file = eg.popup_get_file("Select file")
        layout = setup_layout(load_existing_file(selected_file))
        window = eg.Window(
            "Condition.inp file Generator for Element-based PDD Simulations",
            layout=layout,
        )
    if event == "Cancel":
        window.close()
        eg.popup_warning("condition.inp not generated.")
        break
