# ==============================================================
# incgen.py: Input file generator for incstr code.
#            This script generates inclusions and point defects.
# --------------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# ==============================================================

import TkEasyGUI as eg


class Inclusion:
    def __init__(self):
        self.n = 0
        self.eigen_strains = []
        self.coords = []
        self.sizes = []
        self.orientations = []
        self.parameter = []

    def add(self, eigen_strain, coord, sizes, orientation, parameter):
        self.eigen_strains.append(eigen_strain)
        self.coords.append(coord)
        self.sizes.append(sizes)
        self.orientations.append(orientation)
        self.parameter.append(parameter)
        self.n += 1

    def get_num(self):
        return self.n

    def get_eigen_strain(self, id):
        return self.eigen_strains[id]

    def get_coord(self, id):
        return self.coords[id]

    def get_size(self, id):
        return self.sizes[id]

    def get_orientation(self, id):
        return self.orientations[id]

    def get_parameter(self, id):
        return self.parameter[id]


def open_window():
    spherical_layout = [
        [
            eg.Text("Center position"),
            eg.Text("X"),
            eg.Input("0.0", key="-sphere_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-sphere_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-sphere_z-"),
        ],
        [eg.Text("Radius"), eg.Input("0.0", key="-sphere_radius-")],
        [eg.Text("Eigen strain"), eg.Input("0.0", key="-sphere_eigen_strain-")],
        [eg.HSeparator()],
        [eg.Button("Add", key="-add_spherical_inclusion-")],
    ]
    cylinder_layout = [
        [
            eg.Text("Center position"),
            eg.Text("X"),
            eg.Input("0.0", key="-cylinder_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-cylinder_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-cylinder_z-"),
        ],
        [
            eg.Text("Radius"),
            eg.Input("0.0", key="-cylinder_radius-"),
            eg.Text("Height"),
            eg.Input("0.0", key="-cylinder_height-"),
            eg.Text("1st nnbr distance"),
            eg.Input("2.5e-10", key="-cylinder_1st_nnbr-"),
        ],
        [
            eg.Text("Eigen strain"),
            eg.Input("0.0", key="-cylinder_eigen_strain-"),
            eg.Checkbox("Uniaxial eigen strain", key="-uniaxial_eigenstrain-"),
        ],
        [
            eg.Text("Major axis"),
            eg.Text("X"),
            eg.Input("0.0", key="-cylinder_axis_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-cylinder_axis_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-cylinder_axis_z-"),
        ],
        [eg.HSeparator()],
        [eg.Button("Add", key="-add_cylinderical_inclusion-")],
    ]
    cuboidal_layout = [
        [
            eg.Text("Center position"),
            eg.Text("X"),
            eg.Input("0.0", key="-cuboid_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-cuboid_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-cuboid_z-"),
        ],
        [
            eg.Text("Side length"),
            eg.Text("X"),
            eg.Input("0.0", key="-cuboid_side_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-cuboid_side_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-cuboid_side_z-"),
        ],
        [eg.Text("Eigen strain"), eg.Input("0.0", key="-cuboid_eigen_strain-")],
        [eg.Text("Orientation")],
        [
            eg.Text("  XX"),
            eg.Input("0.0", key="-cuboid_xx-"),
            eg.Text("XY"),
            eg.Input("0.0", key="-cuboid_xy-"),
            eg.Text("XZ"),
            eg.Input("0.0", key="-cuboid_xz-"),
        ],
        [
            eg.Text("  YX"),
            eg.Input("0.0", key="-cuboid_yx-"),
            eg.Text("YY"),
            eg.Input("0.0", key="-cuboid_yy-"),
            eg.Text("YZ"),
            eg.Input("0.0", key="-cuboid_yz-"),
        ],
        [
            eg.Text("  ZX"),
            eg.Input("0.0", key="-cuboid_zx-"),
            eg.Text("ZY"),
            eg.Input("0.0", key="-cuboid_zy-"),
            eg.Text("ZZ"),
            eg.Input("0.0", key="-cuboid_zz-"),
        ],
        [eg.HSeparator()],
        [eg.Button("Add", "-add_cuboidal_inclusion-")],
    ]
    truncated_spherical_layout = [
        [
            eg.Text("Center position"),
            eg.Text("X"),
            eg.Input("0.0", key="-truncated_sphere_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-truncated_sphere_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-truncated_sphere_z-"),
        ],
        [
            eg.Text("Radius"),
            eg.Input("0.0", key="-truncated_sphere_radius-"),
            eg.Text("z1"),
            eg.Input("0.0", key="-truncated_sphere_z1-"),
            eg.Text("z2"),
            eg.Input("0.0", key="-truncated_sphere_z2-"),
        ],
        [
            eg.Text("Eigen strain"),
            eg.Input("0.0", key="-truncated_sphere_eigen_strain-"),
        ],
        [
            eg.Text("Major axis"),
            eg.Text("X"),
            eg.Input("0.0", key="-truncated_sphere_axis_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-truncated_sphere_axis_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-truncated_sphere_axis_z-"),
        ],
        [eg.HSeparator()],
        [eg.Button("Add", "-add_truncated_spherical_inclusion-")],
    ]
    elastic_dipole_layout = [
        [
            eg.Text("Center position"),
            eg.Text("X"),
            eg.Input("0.0", key="-elastic_dipole_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-elastic_dipole_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-elastic_dipole_z-"),
        ],
        [
            eg.Text("Core radius"),
            eg.Input("1.0e-10", key="-elastic_dipole_core_radius-"),
        ],
        [eg.Text("Elastic dipole tensor")],
        [
            eg.Text("  11"),
            eg.Input("0.0", key="-elastic_dipole_11-"),
            eg.Text("12"),
            eg.Input("0.0", key="-elastic_dipole_12-"),
            eg.Text("13"),
            eg.Input("0.0", key="-elastic_dipole_13-"),
        ],
        [
            eg.Text("  21"),
            eg.Input("0.0", key="-elastic_dipole_21-"),
            eg.Text("22"),
            eg.Input("0.0", key="-elastic_dipole_22-"),
            eg.Text("23"),
            eg.Input("0.0", key="-elastic_dipole_23-"),
        ],
        [
            eg.Text("  31"),
            eg.Input("0.0", key="-elastic_dipole_31-"),
            eg.Text("32"),
            eg.Input("0.0", key="-elastic_dipole_32-"),
            eg.Text("33"),
            eg.Input("0.0", key="-elastic_dipole_33-"),
        ],
        [eg.HSeparator()],
        [eg.Button("Add", key="-add_elastic_dipole-")],
    ]
    dilatation_center_layout = [
        [
            eg.Text("Center position"),
            eg.Text("X"),
            eg.Input("0.0", key="-dilatation_center_x-"),
            eg.Text("Y"),
            eg.Input("0.0", key="-dilatation_center_y-"),
            eg.Text("Z"),
            eg.Input("0.0", key="-dilatation_center_z-"),
        ],
        [
            eg.Text("Dilatation strength"),
            eg.Input("0.0", key="-dilatation_center_strength-"),
            eg.Text("Lattice constant"),
            eg.Input("0.0", key="-dilatation_center_a0-"),
        ],
        [eg.HSeparator()],
        [eg.Button("Add", key="-add_dilatational_center-")],
    ]
    layout = [
        [
            eg.TabGroup(
                [
                    eg.Tab(
                        "Spherical",
                        spherical_layout,
                        expand_x=True,
                        expand_y=True,
                    ),
                    eg.Tab(
                        "Cylinderical",
                        cylinder_layout,
                        expand_x=True,
                        expand_y=True,
                    ),
                    eg.Tab(
                        "Cuboidal",
                        cuboidal_layout,
                        expand_x=True,
                        expand_y=True,
                    ),
                    eg.Tab(
                        "Truncated spherical",
                        truncated_spherical_layout,
                        expand_x=True,
                        expand_y=True,
                    ),
                    eg.Tab(
                        "Elastic dipole",
                        elastic_dipole_layout,
                        expand_x=True,
                        expand_y=True,
                    ),
                    eg.Tab(
                        "Dilatation center",
                        dilatation_center_layout,
                        expand_x=True,
                        expand_y=True,
                    ),
                ],
            )
        ],
        [eg.HSeparator()],
        [
            eg.Multiline(
                "0. No inclusion.",
                key="-log-",
                size=(40, 5),
                font=("Arial", 10),
                expand_x=True,
            )
        ],
        [eg.Button("Save"), eg.Button("Cancel")],
    ]
    return eg.Window(
        "Inclusion Generator for incstrgen",
        layout=layout,
        resizable=True,
        size=(900, 500),
    )


def save_inclusions(
    selected_folder,
    sphere,
    cylinder,
    cuboid,
    truncated_sphere,
    elastic_dipole,
    dilatation_center,
):
    with open(f"{selected_folder}/inclusion.inp", "w") as f:
        n = (
            sphere.get_num()
            + cylinder.get_num()
            + cuboid.get_num()
            + truncated_sphere.get_num()
            + elastic_dipole.get_num()
            + dilatation_center.get_num()
        )
        f.write(f"{n}\n")
        for i in range(sphere.get_num()):
            f.write("sphere\n")
            f.write(f"{sphere.get_eigen_strain(i):e}\n")
            f.write(f"{sphere.get_size(i):e}\n")
            coord = sphere.get_coord(i)
            for x in coord:
                f.write(f"{x:e} ")
            f.write("\n")
        for i in range(cylinder.get_num()):
            if cylinder.get_parameter(i) == 1:
                f.write("axial_cylinder\n")
            else:
                f.write("cylinder\n")
            f.write(f"{cylinder.get_eigen_strain(i):e}\n")
            sizes = cylinder.get_size(i)
            for s in sizes:
                f.write(f"{s:e} ")
            f.write("\n")
            coord = cylinder.get_coord(i)
            for x in coord:
                f.write(f"{x:e} ")
            f.write("\n")
            orientation = cylinder.get_orientation(i)
            for x in orientation:
                f.write(f"{x:e} ")
            f.write("\n")
        for i in range(cuboid.get_num()):
            f.write("cuboid\n")
            f.write(f"{cuboid.get_eigen_strain(i):e}\n")
            sizes = cuboid.get_size(i)
            for s in sizes:
                f.write(f"{s:e} ")
            f.write("\n")
            coord = cuboid.get_coord(i)
            for x in coord:
                f.write(f"{x:e} ")
            f.write("\n")
            orientation = cuboid.get_orientation(i)
            for vector in orientation:
                for x in vector:
                    f.write(f"{x:e} ")
                f.write("\n")
        for i in range(truncated_sphere.get_num()):
            f.write("truncated_sphere\n")
            f.write(f"{truncated_sphere.get_eigen_strain(i):e}\n")
            sizes = truncated_sphere.get_size(i)
            for s in sizes:
                f.write(f"{s:e} ")
            f.write("\n")
            coord = truncated_sphere.get_coord(i)
            for x in coord:
                f.write(f"{x:e} ")
            f.write("\n")
            orientation = truncated_sphere.get_orientation(i)
            for x in orientation:
                f.write(f"{x:e} ")
            f.write("\n")
        for i in range(elastic_dipole.get_num()):
            f.write("elastic_dipole\n")
            f.write("0.0 ")
            sizes = elastic_dipole.get_size(i)
            f.write(f"{float(sizes[0]):e}\n")
            coord = elastic_dipole.get_coord(i)
            for x in coord:
                f.write(f"{x:e} ")
            f.write("\n")
            parameter = elastic_dipole.get_parameter(i)
            for vector in parameter:
                for x in vector:
                    f.write(f"{x:e} ")
                f.write("\n")
        for i in range(dilatation_center.get_num()):
            f.write("dilatation_center\n")
            parameter = dilatation_center.get_parameter(i)
            f.write(f"{float(parameter[0]):e}\n")
            f.write(f"{float(parameter[1]):e}\n")
            coord = dilatation_center.get_coord(i)
            for x in coord:
                f.write(f"{x:e} ")
            f.write("\n")


def main():
    sphere = Inclusion()
    cylinder = Inclusion()
    cuboid = Inclusion()
    truncated_sphere = Inclusion()
    elastic_dipole = Inclusion()
    dilatation_center = Inclusion()

    window = open_window()

    log = "0. No inclusion."
    log_count = 1
    while window.is_alive():
        event, values = window.read()
        if event == "-add_spherical_inclusion-":
            eigen_strain = float(values["-sphere_eigen_strain-"])
            coord = [
                float(values["-sphere_x-"]),
                float(values["-sphere_y-"]),
                float(values["-sphere_z-"]),
            ]
            radius = float(values["-sphere_radius-"])
            orientation = []
            parameter = 0
            sphere.add(eigen_strain, coord, radius, orientation, parameter)
            log = (
                f"{log_count}. Spherical inclusion is added: {sphere.get_num()}.\n"
                + log
            )
            log_count += 1
            window["-log-"].update(log)
        if event == "-add_cylinderical_inclusion-":
            eigen_strain = float(values["-cylinder_eigen_strain-"])
            coord = [
                float(values["-cylinder_x-"]),
                float(values["-cylinder_y-"]),
                float(values["-cylinder_z-"]),
            ]
            sizes = [
                float(values["-cylinder_radius-"]),
                float(values["-cylinder_height-"]),
                float(values["-cylinder_1st_nnbr-"]),
            ]
            orientation = [
                float(values["-cylinder_axis_x-"]),
                float(values["-cylinder_axis_y-"]),
                float(values["-cylinder_axis_z-"]),
            ]
            parameter = int(values["-uniaxial_eigenstrain-"])
            cylinder.add(eigen_strain, coord, sizes, orientation, parameter)
            if parameter == 1:
                log = (
                    f"{log_count}. Axial cylindrical inclusion is added: {cylinder.get_num()}.\n"
                    + log
                )
            else:
                log = (
                    f"{log_count}. Cylindrical inclusion is added: {cylinder.get_num()}.\n"
                    + log
                )
            log_count += 1
            window["-log-"].update(log)
        if event == "-add_cuboidal_inclusion-":
            eigen_strain = float(values["-cuboid_eigen_strain-"])
            coord = [
                float(values["-cuboid_x-"]),
                float(values["-cuboid_y-"]),
                float(values["-cuboid_z-"]),
            ]
            sizes = [
                float(values["-cuboid_side_x-"]),
                float(values["-cuboid_side_y-"]),
                float(values["-cuboid_side_z-"]),
            ]
            orientation = [
                [
                    float(values["-cuboid_xx-"]),
                    float(values["-cuboid_xy-"]),
                    float(values["-cuboid_xz-"]),
                ],
                [
                    float(values["-cuboid_yx-"]),
                    float(values["-cuboid_yy-"]),
                    float(values["-cuboid_yz-"]),
                ],
                [
                    float(values["-cuboid_zx-"]),
                    float(values["-cuboid_zy-"]),
                    float(values["-cuboid_zz-"]),
                ],
            ]
            parameter = 0
            cuboid.add(eigen_strain, coord, sizes, orientation, parameter)
            log = (
                f"{log_count}. Cuboidal inclusion is added: {cuboid.get_num()}.\n"
                + log
            )
            log_count += 1
            window["-log-"].update(log)
        if event == "-add_truncated_spherical_inclusion-":
            eigen_strain = float(values["-truncated_sphere_eigen_strain-"])
            coord = [
                float(values["-truncated_sphere_x-"]),
                float(values["-truncated_sphere_y-"]),
                float(values["-truncated_sphere_z-"]),
            ]
            sizes = [
                float(values["-truncated_sphere_radius-"]),
                float(values["-truncated_sphere_z1-"]),
                float(values["-truncated_sphere_z2-"]),
            ]
            orientation = [
                float(values["-truncated_sphere_axis_x-"]),
                float(values["-truncated_sphere_axis_y-"]),
                float(values["-truncated_sphere_axis_z-"]),
            ]
            parameter = 0
            truncated_sphere.add(eigen_strain, coord, sizes, orientation, parameter)
            log = (
                f"{log_count}. Truncated spherical inclusion is added: {truncated_sphere.get_num()}.\n"
                + log
            )
            log_count += 1
            window["-log-"].update(log)
        if event == "-add_elastic_dipole-":
            coord = [                
                float(values["-elastic_dipole_x-"]),
                float(values["-elastic_dipole_y-"]),
                float(values["-elastic_dipole_z-"]),
            ]
            sizes = [
                float(values["-elastic_dipole_core_radius-"])
            ]
            parameter = [
                [
                    float(values["-elastic_dipole_11-"]),
                    float(values["-elastic_dipole_12-"]),
                    float(values["-elastic_dipole_13-"]),
                ],
                [
                    float(values["-elastic_dipole_21-"]),
                    float(values["-elastic_dipole_22-"]),
                    float(values["-elastic_dipole_23-"]),
                ],
                [
                    float(values["-elastic_dipole_31-"]),
                    float(values["-elastic_dipole_32-"]),
                    float(values["-elastic_dipole_33-"]),
                ],
            ]
            eigen_strain = 0
            orientation = 0
            elastic_dipole.add(eigen_strain, coord, sizes, orientation, parameter)
            log = (
                f"{log_count}. Elastic dipole point defect is added: {elastic_dipole.get_num()}.\n"
                + log
            )
            log_count += 1
            window["-log-"].update(log)
        if event == "-add_dilatational_center-":
            coord = [
                float(values["-dilatation_center_x-"]),
                float(values["-dilatation_center_y-"]),
                float(values["-dilatation_center_z-"]),
            ]
            parameter = [
                float(values["-dilatation_center_strength-"]),
                float(values["-dilatation_center_a0-"]),
            ]
            eigen_strain = 0
            sizes = 0
            orientation = 0
            dilatation_center.add(
                eigen_strain, coord, sizes, orientation, parameter
            )
            log = (
                f"{log_count}. Dilatation center point defect is added: {dilatation_center.get_num()}.\n"
                + log
            )
            log_count += 1
            window["-log-"].update(log)
        if event == "Save":
            window.close()
            selected_folder = eg.popup_get_folder("Select folder")
            save_inclusions(
                selected_folder,
                sphere,
                cylinder,
                cuboid,
                truncated_sphere,
                elastic_dipole,
                dilatation_center,
            )
            eg.popup("Inclusions are saved (inclusion.inp).")
            break
        if event == "Cancel":
            window.close()
            eg.popup_warning("Inclusions not saved.")
            break


if __name__ == "__main__":
    main()
