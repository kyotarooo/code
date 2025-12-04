import TkEasyGUI as eg


def write_input(values, selected_folder):
    with open(f"{selected_folder}/section.inp", "w") as f:
        f.write(
            f"{float(values["-volume_x-"]):e} {float(values["-volume_y-"]):e} {float(values["-volume_z-"]):e}\n"
        )
        f.write(f"{float(values["-X-"]):e} {float(values["-Y-"]):e} {float(values["-Z-"]):e}\n")
        f.write(f"{float(values["-NX-"]):e} {float(values["-NY-"]):e} {float(values["-NZ-"]):e}\n")
        f.write(f"{float(values["-FX-"]):e} {float(values["-FY-"]):e} {float(values["-FZ-"]):e}\n")
        f.write(f"{float(values["-PX-"]):e} {float(values["-PY-"]):e} {float(values["-PZ-"]):e}\n")
        f.write(f"{int(values["-CX-"])} {int(values["-CY-"])} {int(values["-CZ-"])}")


def open_window():
    layout = [
        [
            eg.Frame(
                "Volume information",
                [
                    [eg.Text("Volume size")],
                    [
                        eg.Text("X"),
                        eg.Input("1.0", key="-volume_x-"),
                        eg.Text("Y"),
                        eg.Input("1.0", key="-volume_y-"),
                        eg.Text("Z"),
                        eg.Input("1.0", key="-volume_z-"),
                    ],
                    [eg.Text("Number of grid-cells (must be more than 1)")],
                    [
                        eg.Text("NX"),
                        eg.Input("1", key="-CX-"),
                        eg.Text("NY"),
                        eg.Input("1", key="-CY-"),
                        eg.Text("NZ"),
                        eg.Input("1", key="-CZ-"),
                    ],
                ],
                expand_x=True,
            )
        ],
        [
            eg.Frame(
                "Section information",
                [
                    [eg.Text("Position of a point on the section")],
                    [
                        eg.Text("X"),
                        eg.Input("0.0", key="-X-"),
                        eg.Text("Y"),
                        eg.Input("0.0", key="-Y-"),
                        eg.Text("Z"),
                        eg.Input("0.0", key="-Z-"),
                    ],
                    [eg.Text("Normal vector to the section")],
                    [
                        eg.Text("X"),
                        eg.Input("1.0", key="-NX-"),
                        eg.Text("Y"),
                        eg.Input("0.0", key="-NY-"),
                        eg.Text("Z"),
                        eg.Input("0.0", key="-NZ-"),
                    ],
                ],
                expand_x=True,
            )
        ],
        [
            eg.Frame(
                "Stress information",
                [
                    [eg.Text("Force direction")],
                    [
                        eg.Text("X"),
                        eg.Input("0.0", key="-FX-"),
                        eg.Text("Y"),
                        eg.Input("0.0", key="-FY-"),
                        eg.Text("Z"),
                        eg.Input("0.0", key="-FZ-"),
                    ],
                    [eg.Text("Plane direction")],
                    [
                        eg.Text("X"),
                        eg.Input("0.0", key="-PX-"),
                        eg.Text("Y"),
                        eg.Input("0.0", key="-PY-"),
                        eg.Text("Z"),
                        eg.Input("0.0", key="-PZ-"),
                    ],
                ],
                expand_x=True,
            )
        ],
        [eg.Button("Save"), eg.Button("Cancel")],
    ]
    return eg.Window("Input file generator for section code", layout=layout)


def main():
    window = open_window()
    while window.is_alive():
        event, values = window.read()

        if event == "Save":
            window.close()
            selected_folder = eg.popup_get_folder("Select folder")
            if selected_folder:
                write_input(values, selected_folder)
            else:
                eg.popup_warning("Input file not saved.")
            break
        if event == "Cancel":
            window.close()
            eg.popup_warning("Input file not saved.")
            break


if __name__ == "__main__":
    main()
