# =================================================================
# inclusion_vtk.py: vtk file generator from input data for instrgen
# =================================================================
# Usage: python3 vtkgen.py <number_of_vertices>
# [number_of_vertices]
# 1:  Spherical inclusions
# 4:  Cuboidal inclusions
# >4: Cylinderical inclusions
# -----------------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# =================================================================
import sys
import math
import numpy as np


class Template:
    def __init__(self, n):
        self.coords = []

        """Generate vertex coordinates for template based on number of vertices n."""
        if n == 1:
            self.coords.append([0.0, 0.0, 0.0])
        elif n == 4:
            self.coords.append([0.5, 0.5, 0.0])
            self.coords.append([-0.5, 0.5, 0.0])
            self.coords.append([-0.5, -0.5, 0.0])
            self.coords.append([0.5, -0.5, 0.0])
        else:
            dangle = 2.0 * math.pi / n
            angle = 0.5 * dangle
            self.coords = [
                [
                    0.5 * math.cos(angle + i * dangle),
                    0.5 * math.sin(angle + i * dangle),
                    0.0,
                ]
                for i in range(n)
            ]
        self.n = n

    def get_num(self):
        return self.n

    def get_coords(self):
        return self.coords


class Orientation:
    def __init__(self):
        self.tensor = []

    def read_data(self, filename="material.inp"):
        """Read orientation data and normalize each vector."""
        try:
            with open(filename, "r") as f:
                data_lines = f.readlines()
            for i in range(3):
                vector = np.array(list(map(float, data_lines[i].split())))
                vector /= np.linalg.norm(vector)
                self.tensor.append(vector)
        except FileNotFoundError:
            print(f"Error: '{filename}' not found.")
            sys.exit(1)
        except ValueError:
            print("Error: Invalid format in orientation file.")
            sys.exit(1)

    def get_num(self):
        return self.n

    def get_tensor(self):
        return self.tensor


class Inclusion:
    def __init__(self):
        self.n = 0
        self.shapes = []
        self.scales = []
        self.coords = []
        self.tensors = []

    def normalize_vector(self, vector):
        """Normalize a 3D vector."""
        norm = np.linalg.norm(vector)
        if norm == 0:
            raise ValueError("Zero-length vector cannot be normalized.")
        return vector / norm

    def parse_tensor(self, lines):
        """Parse and normalize a 3x3 tensor from lines."""
        tensor = np.array([np.fromstring(line, sep=" ") for line in lines])
        normalized_tensor = np.array([self.normalize_vector(row) for row in tensor])
        return normalized_tensor

    def read_data(self, filename="inclusion.inp"):
        """Read inclusion data, including coordinates, scales, and types."""
        try:
            with open(filename, "r") as f:
                data_lines = f.readlines()

            self.n = int(data_lines[0].strip())
            index = 1

            for _ in range(self.n):
                shape = data_lines[index].strip()
                self.shapes.append(shape)

                if shape == "sphere":
                    radius = data_lines[index + 2].split()[0]
                    self.scales.append(radius)
                    coord = np.fromstring(data_lines[index + 3], sep=" ")
                    self.coords.append(coord[:3])
                    index += 4

                elif shape == "cylinder" or shape == "axial_cylinder":
                    radius, height = map(float, data_lines[index + 2].split())
                    self.scales.append([2.0 * radius, 2.0 * radius, height])
                    coord = np.fromstring(data_lines[index + 3], sep=" ")
                    self.coords.append(coord[:3])
                    axis = np.fromstring(data_lines[index + 4], sep=" ")
                    normalized_axis = self.normalize_vector(axis[:3])
                    tensor = np.array(
                        [
                            normalized_axis,
                            [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0],
                        ]
                    )
                    self.tensors.append(tensor)
                    index += 5

                elif shape == "cuboid":
                    scale = np.fromstring(data_lines[index + 2], sep=" ")
                    self.scales.append(scale[:3])
                    coord = np.fromstring(data_lines[index + 3], sep=" ")
                    self.coords.append(coord[:3])
                    tensor = self.parse_tensor(data_lines[index + 4 : index + 7])
                    self.tensors.append(tensor)
                    index += 7

                elif shape == "elastic_dipole":
                    coord = np.fromstring(data_lines[index + 2], sep=" ")
                    self.coords.append(coord[:3])
                    index += 6

                elif shape == "dilatation_center":
                    a0 = data_lines[index + 2].split()[0]
                    self.scales.append(a0)
                    coord = np.fromstring(data_lines[index + 3], sep=" ")
                    self.coords.append(coord[:3])
                    index += 4

                else:
                    raise ValueError(
                        f"Unknown shape type '{shape}' at line {index + 1}."
                    )

        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.")
            sys.exit(1)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)

    def get_num(self):
        return self.n

    def get_coord(self, id):
        return self.coords[id]

    def get_scale(self, id):
        return self.scales[id]

    def get_shape(self, id):
        return self.shapes[id]

    def get_tensor(self, id):
        return self.tensors[id]

    def set_rotation_tensor(self, orientation):
        rotation_tensor = orientation.get_tensor()
        for i in range(self.n):
            if self.shapes[i] == "cylinder" or self.shapes[i] == "axial_cylinder":
                # Rotate the axis of the cylinder
                rotated_axis = [
                    sum(rotation_tensor[j][k] * self.tensors[i][0][k] for k in range(3))
                    for j in range(3)
                ]

                # Generate a random axis and ensure it's a valid vector
                random_axis = np.array([0.42142421, 0.3423553, 0.532155])
                random_axis /= np.linalg.norm(random_axis)

                # Compute orthogonal vectors using cross products
                vector = np.cross(random_axis, rotated_axis)
                vector /= np.linalg.norm(vector)

                # Construct the rotated tensor
                rotated_tensor = [
                    vector,
                    np.cross(rotated_axis, vector),
                    rotated_axis,
                ]

            elif self.shapes[i] == "cuboid":
                # Get the inclusion tensor
                inclusion_tensor = self.get_tensor(i)

                # Compute the rotated tensor (matrix multiplication)
                rotated_tensor = np.dot(rotation_tensor, inclusion_tensor).tolist()

            # Update the tensor for the current inclusion
            self.tensors[i] = rotated_tensor


def scale_vertices(n_vertices, vertex_coords, scale_factor):
    """Scale vertices based on inclusion's scale factors."""
    z_shift = 0.5 * scale_factor[2]
    return [
        [scale_factor[0] * x, scale_factor[1] * y, z]
        for z in [-z_shift, z_shift]
        for x, y, _ in vertex_coords
    ]


def transform_vertices(n_vertices, vertex_coords, tensor):
    """Transform vertices using the orientation tensor for a given inclusion type."""
    return [
        [sum(tensor[k][j] * vertex_coords[i][k] for k in range(3)) for j in range(3)]
        for i in range(n_vertices)
    ]


def write_vtk(filename, template, inclusion):
    """Write data in VTK format."""
    if template.get_num() == 1:
        n_points = inclusion.get_num()
    else:
        n_points = 2 * template.get_num() * inclusion.get_num()
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Inclusions\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")

        # Write points
        f.write(f"POINTS {n_points} float\n")
        if template.get_num() == 1:
            for i in range(inclusion.get_num()):
                coords = inclusion.get_coord(i)
                for coord in coords:
                    f.write(f"{coord:5e} ")
                f.write("\n")
        else:
            for i in range(inclusion.get_num()):
                scaled_coords = scale_vertices(
                    template.get_num(), template.get_coords(), inclusion.get_scale(i)
                )
                transformed_coords = transform_vertices(
                    2 * template.get_num(),
                    scaled_coords,
                    inclusion.get_tensor(i),
                )
                for coord in transformed_coords:
                    for j in range(3):
                        coord[j] += inclusion.get_coord(i)[j]
                        f.write(f"{coord[j]:5e} ")
                    f.write("\n")

        # Write cells
        if template.get_num() == 1:
            cell_count = inclusion.get_num()
            cell_data_count = 2 * inclusion.get_num()
            f.write(f"CELLS {cell_count} {cell_data_count}\n")
            for i in range(inclusion.get_num()):
                f.write(f"1 {i}\n")
        else:
            cell_count = 3 * inclusion.get_num()
            cell_data_count = inclusion.get_num() * (4 * template.get_num() + 2 + 3)
            f.write(f"CELLS {cell_count} {cell_data_count}\n")
            for i in range(inclusion.get_num()):
                n_vertices_offset = 2 * template.get_num() * i
                f.write(
                    f"{template.get_num()} {' '.join(map(str, range(n_vertices_offset, n_vertices_offset + template.get_num())))}\n"
                )
                f.write(
                    f"{template.get_num()} {' '.join(map(str, range(n_vertices_offset + template.get_num(), n_vertices_offset + 2 * template.get_num())))}\n"
                )
                f.write(f"{2 * template.get_num() + 2} ")
                for j in range(template.get_num()):
                    f.write(
                        f"{n_vertices_offset + j} {n_vertices_offset + template.get_num() + j} "
                    )
                f.write(
                    f"{n_vertices_offset} {n_vertices_offset + template.get_num()}\n"
                )

        # Write cell types
        f.write(f"CELL_TYPES {cell_count}\n")
        if template.get_num() == 1:
            for i in range(inclusion.get_num()):
                f.write("1\n")  # VTK_VERTEX
        else:
            for i in range(inclusion.get_num()):
                f.write("7\n" * 2)  # VTK_POLYGON
                f.write("6\n")  # VTK_TRIANGLE_STRIP

        # Write radius if inclusions are spherical
        if inclusion.get_shape(0) == "sphere":
            f.write(f"POINT_DATA {inclusion.get_num()}\n")
            f.write("SCALARS radius float\n")
            f.write("LOOKUP_TABLE default\n")
            for i in range(inclusion.get_num()):
                f.write(f"{inclusion.scales[i]}\n")


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <number_of_vertices>")
        print("[number_of_vertices]")
        print("1:  Spherical inclusions/Elastic dipole/Dilatation center")
        print("4:  Cuboidal inclusions")
        print(">4: Cylinderical inclusions")
        sys.exit(1)

    try:
        n_vertices = int(sys.argv[1])
    except ValueError:
        print("Error: Number of vertices must be an integer.")
        sys.exit(1)

    template = Template(n_vertices)

    inclusion = Inclusion()
    inclusion.read_data()
    if template.get_num() != 1:
        orientation = Orientation()
        orientation.read_data()
        inclusion.set_rotation_tensor(orientation)

    write_vtk("inclusion.vtk", template, inclusion)


if __name__ == "__main__":
    main()
