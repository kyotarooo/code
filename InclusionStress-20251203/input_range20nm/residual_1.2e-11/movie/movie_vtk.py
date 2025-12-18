import numpy as np
import sys
import subprocess


class Movie:
    def __init__(self):
        self.n = 1
        self.start_coord = np.zeros(3)
        self.end_coord = np.zeros(3)


class Section:
    def __init__(self):
        self.x = np.zeros(3)
        self.n = np.zeros(3)
        self.s = np.zeros(3)
        self.m = np.zeros(3)


class Volume:
    def __init__(self):
        self.section = Section()
        self.size = np.zeros(3)
        self.nx = np.zeros(3, dtype=int)


def read_movie(movie):
    with open(f"{sys.argv[1]}/movie.inp", "r") as f:
        line = f.readline()
        movie.n = int(line.split()[0])
        line = f.readline()
        movie.start_coord = np.array([float(line.split()[i]) for i in range(3)])
        line = f.readline()
        movie.end_coord = np.array([float(line.split()[i]) for i in range(3)])


def read_section(volume):
    with open(f"{sys.argv[1]}/section.inp", "r") as f:
        lines = f.readlines()
        volume.size = np.array([float(lines[0].split()[i]) for i in range(3)])
        volume.section.x = np.array([float(lines[1].split()[i]) for i in range(3)])
        volume.section.n = np.array([float(lines[2].split()[i]) for i in range(3)])
        volume.section.s = np.array([float(lines[3].split()[i]) for i in range(3)])
        volume.section.m = np.array([float(lines[4].split()[i]) for i in range(3)])
        volume.nx = np.array([int(lines[5].split()[i]) for i in range(3)])


def section_position(i, movie, volume):
    dcoord = (movie.end_coord - movie.start_coord) / (movie.n - 1.0)
    volume.section.x = movie.start_coord + i * dcoord


def write_section(volume):
    size = volume.size
    x = volume.section.x
    n = volume.section.n
    s = volume.section.s
    m = volume.section.m
    nx = volume.nx
    with open(f"{sys.argv[1]}/section.inp", "w") as f:
        f.write(f"{size[0]} {size[1]} {size[2]}\n")
        f.write(f"{x[0]} {x[1]} {x[2]}\n")
        f.write(f"{n[0]} {n[1]} {n[2]}\n")
        f.write(f"{s[0]} {s[1]} {s[2]}\n")
        f.write(f"{m[0]} {m[1]} {m[2]}\n")
        f.write(f"{nx[0]} {nx[1]} {nx[2]}\n")


def take_sceen(i):
    command = "./section_vtk " + sys.argv[1]
    subprocess.run(command, shell=True)
    file_name = f"{sys.argv[1]}/section_{format(i + 1, "04")}.vtk"
    command = "mv " + f"{sys.argv[1]}/section.vtk " + file_name
    subprocess.run(command, shell=True)


def main():
    movie = Movie()
    volume = Volume()

    read_movie(movie)
    read_section(volume)

    for i in range(movie.n):
        section_position(i, movie, volume)
        write_section(volume)
        take_sceen(i)


if __name__ == "__main__":
    main()
