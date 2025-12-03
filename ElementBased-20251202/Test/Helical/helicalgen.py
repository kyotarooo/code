import numpy as np

# Parameters
n = 40
volume_length = 20.0e-09
burgers_vector = 2.5e-10
radius = 5.0e-09
straight_length = 2.0e-09

# Geometry setup
x = 0.5 * straight_length
dx = (volume_length - straight_length) / n
y0 = 0.5 * volume_length
z0 = 0.5 * volume_length

dangle = 2.0 * np.pi / n
angle = 0.0

coords = []

# Write coordinates and element information
with open("elem.inp", "w") as f:
    # Node count
    f.write(f"{n + 1}\n")

    # Node coordinates
    for i in range(n + 1):
        y = y0 + radius * np.sin(angle)
        z = z0 + radius * np.cos(angle)

        f.write(f"2 {x:.6e} {y:.6e} {z:.6e}\n")
        coords.append([x, y, z])

        x += dx
        angle += dangle

    # Element count
    f.write(f"{n + 1}\n")

    # Tangent and normal vectors
    x_direction = np.array([1.0, 0.0, 0.0])
    for i in range(n):
        tangent = np.array(coords[i + 1]) - np.array(coords[i])
        a = np.cross(tangent, x_direction)
        normal = a / np.linalg.norm(a)

        f.write(
            f"{i} {i + 1} {burgers_vector:.6e} 0.0 0.0 "
            f"{normal[0]:.6e} {normal[1]:.6e} {normal[2]:.6e}\n"
        )

    # Final line (closing element)
    f.write(f"{n} 0 {burgers_vector:.6e} 0.0 0.0 0.0 1.0 0.0\n")
