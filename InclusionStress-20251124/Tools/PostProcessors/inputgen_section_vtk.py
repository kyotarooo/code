import TkEasyGUI as eg


def write_input(values, selected_folder):
    with open(f"{selected_folder}/section.inp", "w") as f:
        f.write(
            f"{values["-volume_x-"]} {values["-volume_y-"]} {values["-volume_z-"]}\n"
        )
        f.write(f"{values["-X-"]} {values["-Y-"]} {values["-Z-"]}\n")
        f.write(f"{values["-NX-"]} {values["-NY-"]} {values["-NZ-"]}\n")
        f.write(f"{values["-FX-"]} {values["-FY-"]} {values["-FZ-"]}\n")
        f.write(f"{values["-PX-"]} {values["-PY-"]} {values["-PZ-"]}\n")
        f.write(f"{int(values["-CX-"])} {int(values["-CY-"])} {int(values["-CZ-"])}")

def open_window():
    layout = [
        [eg.Text("Volume information", background_color="gray", text_color="white")],
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
            eg.Text("NX"), eg.Input("1", key="-CX-"),
            eg.Text("NY"), eg.Input("1", key="-CY-"),
            eg.Text("NZ"), eg.Input("1", key="-CZ-"),
         ],
        [eg.HSeparator()],
        [eg.Text("Section information", background_color="gray", text_color="white")],
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
        [eg.HSeparator()],
        [eg.Text("Stress information", background_color="gray", text_color="white")],
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
        [eg.HSeparator()],
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
                eg.popup("Input file not saved.")
            break
        if event == "Cancel":
            window.close()
            eg.popup("Input file not saved.")
            break


if __name__ == "__main__":
    main()
