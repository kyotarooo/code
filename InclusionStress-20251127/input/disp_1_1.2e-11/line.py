import numpy as np
import sys
import struct


class Grid:
    def __init__(self):
        self.nx = np.zeros(3, dtype=int)
        self.size = np.zeros(3)
        self.stress = None


class Line:
    def __init__(self):
        self.n = 0
        self.x0 = np.zeros(3)
        self.x1 = np.zeros(3)
        self.s = np.zeros(3)
        self.m = np.zeros(3)


class Volume:
    def __initi__(self):
        self.size = np.zeros(3)
        self.nx = np.zeros(3, dtype=int)
        self.grid = Grid()
        self.line = Line()


def unit_vector(a):
    r = np.linalg.norm(a)
    a /= r


def read_volume(volume):
    with open(f"{sys.argv[1]}/line.inp", "r") as f:
        lines = f.readlines()
        volume.size = np.array([float(lines[0].split()[i]) for i in range(3)])
        volume.n = int(lines[1].split()[0])
        volume.line.x0 = np.array([float(lines[2].split()[i]) for i in range(3)])
        volume.line.x1 = np.array([float(lines[3].split()[i]) for i in range(3)])
        volume.line.s = np.array([float(lines[4].split()[i]) for i in range(3)])
        volume.line.m = np.array([float(lines[5].split()[i]) for i in range(3)])
        volume.nx = np.array([int(lines[6].split()[i]) for i in range(3)])

        unit_vector(volume.line.s)
        unit_vector(volume.line.m)


def stress_component(line, local_stress):
    index = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
    v = 0.0
    for i in range(3):
        for j in range(3):
            v += local_stress[index[i][j]] * line.s[i] * line.m[j]
    return v


def read_stress(volume):
    with open(f"{sys.argv[1]}/grid_stress.inp", "rb") as f:
        volume.grid.nx = np.array(struct.unpack("3i", f.read(12)))
        volume.grid.size = volume.size / volume.nx / (volume.grid.nx - 1.0)

        stress_shape = tuple(volume.grid.nx)
        stress = np.zeros(stress_shape)

        for i in range(volume.grid.nx[0]):
            for j in range(volume.grid.nx[1]):
                for k in range(volume.grid.nx[2]):
                    active = struct.unpack("i", f.read(4))[0]
                    local_stress = struct.unpack("9d", f.read(72))
                    stress[i, j, k] = stress_component(volume.line, local_stress)

def offset(volume, x):
    for i in range(3):
        size = volume.size[i] / volume.nx[i]
        id = int(math.floor(x[i] / size))
        id = max(0, min(id, volume.nx[i] - 1))
        x[i] -= id * size


def shape_function(u):
    ux, uy, uz = u
    return np.array(
        [
            0.125 * (1 - ux) * (1 - uy) * (1 - uz),
            0.125 * (1 + ux) * (1 - uy) * (1 - uz),
            0.125 * (1 + ux) * (1 + uy) * (1 - uz),
            0.125 * (1 - ux) * (1 + uy) * (1 - uz),
            0.125 * (1 - ux) * (1 - uy) * (1 + uz),
            0.125 * (1 + ux) * (1 - uy) * (1 + uz),
            0.125 * (1 + ux) * (1 + uy) * (1 + uz),
            0.125 * (1 - ux) * (1 + uy) * (1 + uz),
        ]
    )

def interpolate_stress(volume, x):
    grid = volume.grid
    grid_id = np.zeros(3, dtype=int)
    u = np.zeros(3)
    local_stress = 0.0

    offset(volume, x)

    for i in range(3):
        grid_id[i] = int(x[i] / grid.size[i])
        grid_id[i] = max(0, min(grid_id[i], grid.nx[i] - 2))
        u[i] = 2.0 * (x[i] / grid.size[i] - grid_id[i]) - 1.0
        u[i] = max(-1.0, min(u[i], 1.0))

    shape = shape_function(u)

    grid_index = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    )

    for i in range(8):
        idx = tuple(grid_id + grid_index[i])
        local_stress += shape[i] * grid.stress[idx]
    return local_stress

def write_stress_on_line(volume):
    dx = (volume.line.x1 - volume.line.x0) / (volume.n - 1.0)
    with open(f"{sys.argv[1]}/line_stress.txt", "w") as f:
        for i in range(volume.n):
            x = volume.line.x0 + i * dx
            stress = interpolate_stress(volume, x)

            for xi in x:
                f.write(f"{xi} ")
            f.write(f"{stress}\n")
          
def main():
    volume = Volume()

    read_volume(volume)
    read_stress(volume)

    write_stress_on_line(volume)

if __name__ == "__main__":
    main()
