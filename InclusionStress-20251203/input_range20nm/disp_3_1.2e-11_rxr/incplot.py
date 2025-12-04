import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

SPHERE_RESOLUTION = 8
CYLINDER_RESOLUTION = 8

fig = plt.figure("Inclusion plot")
ax = fig.add_subplot(projection="3d")


class InclusionPlot:
    def __init__(self):
        n = SPHERE_RESOLUTION
        u = np.linspace(0.0, 2.0 * np.pi, n)
        v = np.linspace(0.0, np.pi, n)
        self.templete_sphere_x = 0.5 * np.outer(np.cos(u), np.sin(v))
        self.templete_sphere_y = 0.5 * np.outer(np.sin(u), np.sin(v))
        self.templete_sphere_z = 0.5 * np.outer(np.ones(np.size(u)), np.cos(v))

        n = CYLINDER_RESOLUTION
        rho = np.linspace(0.0, 0.5, 2)
        phi = np.linspace(0.0, 2.0 * np.pi, n)
        self.templete_cylinder_top_x = np.outer(rho, np.cos(phi))
        self.templete_cylinder_top_y = np.outer(rho, np.sin(phi))
        self.templete_cylinder_top_z = np.outer(0.5 * np.ones(n), np.ones(n))

        self.templete_cylinder_bottom_x = np.outer(rho, np.cos(phi))
        self.templete_cylinder_bottom_y = np.outer(rho, np.sin(phi))
        self.templete_cylinder_bottom_z = np.outer(-0.5 * np.ones(n), np.ones(n))

        self.templete_cylinder_side_x = np.outer(0.5 * np.ones(2), np.cos(phi))
        self.templete_cylinder_side_y = np.outer(0.5 * np.ones(2), np.sin(phi))
        z = np.linspace(-0.5, 0.5, 2)
        self.templete_cylinder_side_z = np.outer(z, np.ones(n))

    def init_sphere(self):
        self.sphere_x = self.templete_sphere_x.copy()
        self.sphere_y = self.templete_sphere_y.copy()
        self.sphere_z = self.templete_sphere_z.copy()

    def scale_sphere(self, scale):
        self.sphere_x *= scale
        self.sphere_y *= scale
        self.sphere_z *= scale

    def translate_sphere(self, translation):
        self.sphere_x += translation[0]
        self.sphere_y += translation[1]
        self.sphere_z += translation[2]

    def plot_sphere(self):
        ax.plot_surface(
            self.sphere_x,
            self.sphere_y,
            self.sphere_z,
            color="lightgray",
            antialiased=False,
        )

    def init_cylinder(self):
        self.cylinder_top_x = self.templete_cylinder_top_x.copy()
        self.cylinder_top_y = self.templete_cylinder_top_y.copy()
        self.cylinder_top_z = self.templete_cylinder_top_z.copy()
        self.cylinder_bottom_x = self.templete_cylinder_bottom_x.copy()
        self.cylinder_bottom_y = self.templete_cylinder_bottom_y.copy()
        self.cylinder_bottom_z = self.templete_cylinder_bottom_z.copy()
        self.cylinder_side_x = self.templete_cylinder_side_x.copy()
        self.cylinder_side_y = self.templete_cylinder_side_y.copy()
        self.cylinder_side_z = self.templete_cylinder_side_z.copy()

    def scale_cylinder(self, diameter_scale, height_scale):
        self.cylinder_top_x *= diameter_scale
        self.cylinder_top_y *= diameter_scale
        self.cylinder_bottom_x *= diameter_scale
        self.cylinder_bottom_y *= diameter_scale
        self.cylinder_side_x *= diameter_scale
        self.cylinder_side_y *= diameter_scale

        self.cylinder_top_z *= height_scale
        self.cylinder_bottom_z *= height_scale
        self.cylinder_side_z *= height_scale

    def rotate_coordinate(self, tensor, coord):
        return [sum(tensor[j][i] * coord[j] for j in range(3)) for i in range(3)]

    def rotate_cylinder(self, material, direction):
        rotation_tensor = material.get_tensor()
        rotated_axis = [
            sum(rotation_tensor[i][j] * direction[j] for k in range(3))
            for i in range(3)
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

        n = CYLINDER_RESOLUTION
        for i in range(2):
            for j in range(n):
                coord = np.array(
                    [
                        self.cylinder_top_x[i][j],
                        self.cylinder_top_y[i][j],
                        self.cylinder_top_z[i][j],
                    ]
                )
                rotated_coord = self.rotate_coordinate(rotated_tensor, coord)
                self.cylinder_top_x[i][j] = rotated_coord[0]
                self.cylinder_top_y[i][j] = rotated_coord[1]
                self.cylinder_top_z[i][j] = rotated_coord[2]

                coord = np.array(
                    [
                        self.cylinder_bottom_x[i][j],
                        self.cylinder_bottom_y[i][j],
                        self.cylinder_bottom_z[i][j],
                    ]
                )
                rotated_coord = self.rotate_coordinate(rotated_tensor, coord)
                self.cylinder_bottom_x[i][j] = rotated_coord[0]
                self.cylinder_bottom_y[i][j] = rotated_coord[1]
                self.cylinder_bottom_z[i][j] = rotated_coord[2]

                coord = np.array(
                    [
                        self.cylinder_side_x[i][j],
                        self.cylinder_side_y[i][j],
                        self.cylinder_side_z[i][j],
                    ]
                )
                rotated_coord = self.rotate_coordinate(rotated_tensor, coord)
                self.cylinder_side_x[i][j] = rotated_coord[0]
                self.cylinder_side_y[i][j] = rotated_coord[1]
                self.cylinder_side_z[i][j] = rotated_coord[2]

    def translate_cylinder(self, translation):
        self.cylinder_top_x += translation[0]
        self.cylinder_top_y += translation[1]
        self.cylinder_top_z += translation[2]
        self.cylinder_bottom_x += translation[0]
        self.cylinder_bottom_y += translation[1]
        self.cylinder_bottom_z += translation[2]
        self.cylinder_side_x += translation[0]
        self.cylinder_side_y += translation[1]
        self.cylinder_side_z += translation[2]

    def plot_cylinder(self):
        ax.plot_surface(
            self.cylinder_side_x,
            self.cylinder_side_y,
            self.cylinder_side_z,
            color="lightgray",
            antialiased=False,
        )
        ax.plot_surface(
            self.cylinder_bottom_x,
            self.cylinder_bottom_y,
            self.cylinder_bottom_z,
            color="lightgray",
            antialiased=False,
        )
        ax.plot_surface(
            self.cylinder_top_x,
            self.cylinder_top_y,
            self.cylinder_top_z,
            color="lightgray",
            antialiased=False,
        )


class Material:
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
            self.size = np.array(list(map(float, data_lines[4].split())))
        except FileNotFoundError:
            print(f"Error: '{filename}' not found.")
            sys.exit(1)
        except ValueError:
            print("Error: Invalid format in orientation file.")
            sys.exit(1)

    def get_tensor(self):
        return self.tensor

    def get_size(self):
        return self.size


class Inclusion:
    def read_and_plot(self, inclusion_plot, material):
        with open("inclusion.inp", "r") as f:
            data_lines = f.readlines()
        n = int(data_lines[0].strip())
        index = 1
        for i in range(n):
 #           shape = data_lines[index].strip()
            shape = data_lines[index].split()[0]
            if shape == "sphere":
 #               radius = float(data_lines[index + 2].strip())
 #               coord = np.array(list(map(float, data_lines[index + 3].split())))
 #               inclusion_plot.init_sphere()
 #               inclusion_plot.scale_sphere(2.0 * radius)
 #               inclusion_plot.translate_sphere(coord[:3])
 #               inclusion_plot.plot_sphere()
 #               index = index + 4
                radius = float(data_lines[index].split()[2])
                coord = [float(data_lines[index].split()[3]), float(data_lines[index].split()[4]), float(data_lines[index].split()[5])]
                inclusion_plot.init_sphere()
                inclusion_plot.scale_sphere(2.0 * radius)
                inclusion_plot.translate_sphere(coord)
                inclusion_plot.plot_sphere()
                index = index + 1
            elif shape == "cylinder":
                radius = float(data_lines[index + 2].split()[0])
                height = float(data_lines[index + 2].split()[1])
                coord = np.array(list(map(float, data_lines[index + 3].split())))
                direction = np.array(list(map(float, data_lines[index + 4].split())))
                direction /= np.linalg.norm(direction)
                inclusion_plot.init_cylinder()
                inclusion_plot.scale_cylinder(2.0 * radius, height)
                inclusion_plot.rotate_cylinder(material, direction[:3])
                inclusion_plot.translate_cylinder(coord[:3])
                index = index + 5


def set_up_plot(material):
    size = material.get_size()

    dx, dy, dz = ([0.0, size[0]], [0.0, size[1]], [0.0, size[2]])

    ax.set_xlim(0.0, size[0])
    ax.set_ylim(0.0, size[1])
    ax.set_zlim(0.0, size[2])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    aspect_ratios = [size[i] / size[0] for i in range(3)]
    ax.set_box_aspect(aspect_ratios)

    ax.grid(False)

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_major_locator(ticker.MaxNLocator(5))
        axis.set_minor_locator(ticker.MaxNLocator(10))
        axis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

    ax.plot([dx[0], dx[1]], [dy[0], dy[0]], [dz[0], dz[0]], color="black", lw=0.5)
    ax.plot([dx[1], dx[1]], [dy[0], dy[1]], [dz[0], dz[0]], color="black", lw=0.5)
    ax.plot([dx[1], dx[0]], [dy[1], dy[1]], [dz[0], dz[0]], color="black", lw=0.5)
    ax.plot([dx[0], dx[0]], [dy[1], dy[0]], [dz[0], dz[0]], color="black", lw=0.5)

    ax.plot([dx[0], dx[0]], [dy[0], dy[0]], [dz[0], dz[1]], color="black", lw=0.5)
    ax.plot([dx[1], dx[1]], [dy[0], dy[0]], [dz[0], dz[1]], color="black", lw=0.5)
    ax.plot([dx[1], dx[1]], [dy[1], dy[1]], [dz[0], dz[1]], color="black", lw=0.5)
    ax.plot([dx[0], dx[0]], [dy[1], dy[1]], [dz[0], dz[1]], color="black", lw=0.5)

    ax.plot([dx[0], dx[1]], [dy[0], dy[0]], [dz[1], dz[1]], color="black", lw=0.5)
    ax.plot([dx[1], dx[1]], [dy[0], dy[1]], [dz[1], dz[1]], color="black", lw=0.5)
    ax.plot([dx[1], dx[0]], [dy[1], dy[1]], [dz[1], dz[1]], color="black", lw=0.5)
    ax.plot([dx[0], dx[0]], [dy[1], dy[0]], [dz[1], dz[1]], color="black", lw=0.5)


def main():
    inclusion_plot = InclusionPlot()
    material = Material()
    inclusion = Inclusion()

    material.read_data()
    set_up_plot(material)
    inclusion.read_and_plot(inclusion_plot, material)

    plt.show()


if __name__ == "__main__":
    main()
