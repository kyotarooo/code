# =====================================================================
# elemgen.py: elem.inp generator for straight dislocations with GUI
# ---------------------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# =====================================================================
import numpy as np
import TkEasyGUI as eg

ZERO_BURGERS_VECTOR = 1
ZERO_SLIP_PLANE = 2
ZERO_ELEMENT = 3
ZERO_NELEMENTS = 4
NOT_ON_SLIP_PLANE = 5


class Dislocations:
    def __init__(self):
        self.npoints0 = 0
        self.nelements0 = 0
        self.points = []
        self.point_types = []
        self.elements_nodes = []
        self.burgers_vectors = []
        self.slip_planes = []

    def set_parameters(
        self,
        start_node,
        end_node,
        burgers_vector,
        burgers_magnitude,
        slip_plane,
        nelements,
        loop_type,
    ):
        self.start_node = start_node
        self.end_node = end_node
        self.burgers_vector = burgers_vector
        self.burgers_magnitude = burgers_magnitude
        self.slip_plane = slip_plane
        self.nelements = nelements
        self.loop_type = loop_type

    def generate_dislocation(self):
        if self.nelements == 0:
            return ZERO_NELEMENTS

        element = [self.end_node[i] - self.start_node[i] for i in range(3)]
        norm = np.linalg.norm(element)
        if norm == 0:
            return ZERO_ELEMENT
        element = element / norm
        for i in range(self.nelements + 1):
            point = [
                self.start_node[j]
                + (self.end_node[j] - self.start_node[j]) * (i / self.nelements)
                for j in range(3)
            ]
            self.points.append(point)
            self.point_types.append(0)

        norm = np.linalg.norm(self.burgers_vector)
        if norm == 0:
            return ZERO_BURGERS_VECTOR
        self.burgers_vector = self.burgers_vector / norm * self.burgers_magnitude

        for i in range(self.nelements):
            self.elements_nodes.append([self.npoints0 + i, self.npoints0 + i + 1])
            self.burgers_vectors.append(self.burgers_vector)
            self.slip_planes.append(self.slip_plane)

        norm = np.linalg.norm(self.slip_plane)
        if norm == 0:
            return ZERO_SLIP_PLANE
        self.slip_plane = self.slip_plane / norm

        cos_angle = abs(sum([self.slip_plane[i] * element[i] for i in range(3)]))
        if cos_angle > 0.001:
            return NOT_ON_SLIP_PLANE

        if self.loop_type == "close":
            self.elements_nodes.append([self.npoints0 + self.nelements, self.npoints0])
            self.burgers_vectors.append(self.burgers_vector)
            self.slip_planes.append(self.slip_plane)
            self.nelements += 1
            self.npoints0 += self.nelements
            self.nelements0 += self.nelements

        else:
            self.point_types[self.npoints0] = 1
            self.point_types[self.npoints0 + self.nelements] = 1
            self.npoints0 += self.nelements + 1
            self.nelements0 += self.nelements

        return 0

    def write_to_inp(self, selected_folder):
        with open(f"{selected_folder}/elem.inp", "w") as f:
            f.write(f"{self.npoints0}\n")
            for i, point in enumerate(self.points):
                line = (
                    f"{self.point_types[i]} "
                    + " ".join(f"{coord:e}" for coord in point)
                    + "\n"
                )
                f.write(line)

            f.write(f"{self.nelements0}\n")
            for i, element in enumerate(self.elements_nodes):
                line = (
                    " ".join(map(str, element))
                    + " "
                    + " ".join(f"{val:e}" for val in self.burgers_vectors[i])
                    + " "
                    + " ".join(f"{val:e}" for val in self.slip_planes[i])
                    + "\n"
                )
                f.write(line)

    def write_to_vtk(self, selected_folder):
        with open(f"{selected_folder}/elem.vtk", "w") as f:
            f.write("vtk DataFile Version 3.0\n")
            f.write("0.0\n")
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")
            f.write(f"POINTS {self.npoints0} float\n")
            for i, point in enumerate(self.points):
                line = " ".join(f"{coord:e}" for coord in point) + "\n"
                f.write(line)
            f.write(f"CELLS {self.nelements0} {self.nelements0 * 3}\n")
            for i, element in enumerate(self.elements_nodes):
                line = "2 " + " ".join(map(str, element)) + "\n"
                f.write(line)
            f.write(f"CELL_TYPES {self.nelements0}\n")
            for i in range(self.nelements0):
                f.write("3\n")
            f.write(f"POINT_DATA {self.npoints0}\n")
            f.write("SCALARS NodeType float\n")
            f.write("LOOKUP_TABLE default\n")
            for type in self.point_types:
                f.write(f"{type:d}\n")
            f.write(f"CELL_DATA {self.nelements0}\n")
            f.write("VECTORS Burgers float\n")
            for i, burgers_vector in enumerate(self.burgers_vectors):
                line = " ".join(f"{val:e}" for val in burgers_vector) + "\n"
                f.write(line)
            f.write("VECTORS SlipPlane float\n")
            for i, slip_plane in enumerate(self.slip_planes):
                line = " ".join(f"{val:e}" for val in slip_plane) + "\n"
                f.write(line)


def main():
    log_count = 0
    log = f"{log_count}. No dislocation.\n"
    log_count += 1
    layout = [
        [eg.Text("Position of node ID:0 (m)")],
        [
            eg.Text("X"),
            eg.Input("0.0", key="-NODEID0X-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-NODEID0Y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-NODEID0Z-"),
        ],
        [eg.Text("Position of node ID:1 (m)")],
        [
            eg.Text("X"),
            eg.Input("0.0", key="-NODEID1X-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-NODEID1Y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-NODEID1Z-"),
        ],
        [eg.Text("Burgers vector direction")],
        [
            eg.Text("X"),
            eg.Input("0.0", key="-BURGERSX-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-BURGERSY-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-BURGERSZ-"),
        ],
        [eg.Text("Burgers vector magnitude (m)"), eg.Input("2.5e-10", key="-BURGERSM-")],
        [eg.Text("Slip plane normal direction")],
        [
            eg.Text("X"),
            eg.Input("0.0", key="-SLIPPLANEX-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-SLIPPLANEY-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-SLIPPLANEZ-"),
        ],
        [eg.Text("Number of elements"), eg.Input("10", key="-NELEMENTS-")],
        [
            eg.Text("Dislocation type"),
            eg.Radio("Open", default=True, group_id=0, key="-OPEN-"),
            eg.Radio("Close", default=False, group_id=0, key="-CLOSE-"),
        ],
        [eg.HSeparator()],
        [
            eg.Multiline(
                log, key="-log-", size=(40, 5), font=("Arial", 10), expand_x=True
            )
        ],
        [eg.Button("Add"), eg.Button("Save"), eg.Button("Cancel")],
    ]
    window = eg.Window(
        "elem.inp file Generator for Element-based PDD Simulations", layout=layout
    )
    generator = Dislocations()

    nloops = 0
    while window.is_alive():
        event, values = window.read()
        if event == "Add":
            node0 = [
                float(values["-NODEID0X-"]),
                float(values["-NODEID0Y-"]),
                float(values["-NODEID0Z-"]),
            ]
            node1 = [
                float(values["-NODEID1X-"]),
                float(values["-NODEID1Y-"]),
                float(values["-NODEID1Z-"]),
            ]
            burgers_vector = [
                float(values["-BURGERSX-"]),
                float(values["-BURGERSY-"]),
                float(values["-BURGERSZ-"]),
            ]
            burgers_magnitude = float(values["-BURGERSM-"])
            slip_plane = [
                float(values["-SLIPPLANEX-"]),
                float(values["-SLIPPLANEY-"]),
                float(values["-SLIPPLANEZ-"]),
            ]
            nelements = int(values["-NELEMENTS-"])
            loop_type = "open" if values["-OPEN-"] else "close"

            generator.set_parameters(
                node0,
                node1,
                burgers_vector,
                burgers_magnitude,
                slip_plane,
                nelements,
                loop_type,
            )
            result = generator.generate_dislocation()
            if result == ZERO_BURGERS_VECTOR:
                eg.popup("Zero Burgers vector cannot be normalized.")
            elif result == ZERO_SLIP_PLANE:
                eg.popup("Zero slip plane vector cannot be normalized.")
            elif result == ZERO_ELEMENT:
                eg.popup("Dislocation length must be non-zero.")
            elif result == ZERO_NELEMENTS:
                eg.popup("Number of elements must be non-zero.")
            elif result == NOT_ON_SLIP_PLANE:
                eg.popup("Dislocation must be on the slip plane.")
            else:
                nloops += 1
                log = f"{log_count}. Dislocation loop is added: {nloops}.\n" + log
                log_count += 1
                window["-log-"].update(log)
        if event == "Save":
            window.close()

            if nloops == 0:
                eg.popup("No dislocation loop is generated.")
            else:
                selected_folder = eg.popup_get_folder("Select folder")
                generator.write_to_inp(selected_folder)
                generator.write_to_vtk(selected_folder)
                eg.popup("elem.inp file is generated.")
                break
        if event == "Cancel":
            window.close()
            eg.popup("elem.inp file not generated.")
            break


if __name__ == "__main__":
    main()
