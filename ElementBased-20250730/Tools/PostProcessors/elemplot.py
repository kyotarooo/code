# ========================================================================================
# elemplot.py: Simple result visualizer for Element-based Parametric Dislocation Dynamics
# ----------------------------------------------------------------------------------------
# Usage:
# - "Up" key: Plot the next step
# - "Down" key: Plot the previous step
# - "n" key: Plot / Not plot nodal points
# - "l" key: Plot / Not plot lines
# - "a" key: Plot / Not plot points and lines outside the simulation volume
# - "i" key: Save figure as image file
# - "m" key: Save images for movie
# - "q" key: Exit
# ----------------------------------------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# ========================================================================================
import os
import sys
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import vtk
from tqdm import tqdm
import cv2

# Please modify as you like
LINE_WIDTH = 1.0
LINE_COLOR = "black"
POINT_SIZE = 10
POINT_COLOR = "white"
POINT_EDGE_COLOR = "black"
FPS = 30

fig = plt.figure("Element-based Parametric Dislocation Dynamics")
ax = fig.add_subplot(projection="3d")


class ElementPlotter:
    def __init__(
        self,
        nsteps,
        line_width=LINE_WIDTH,
        line_color=LINE_COLOR,
        point_size=POINT_SIZE,
        point_color=POINT_COLOR,
        point_edge_color=POINT_EDGE_COLOR,
        fps=FPS,
    ):
        self.line_width = line_width
        self.line_color = line_color
        self.point_size = point_size
        self.point_color = point_color
        self.point_edge_color = point_edge_color
        self.fps = fps
        self.step = 1
        self.plot_points = True
        self.plot_lines = True
        self.plot_outside = False
        self.nsteps = nsteps

    def read_volume(self):
        try:
            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName("volume.vtk")
            reader.Update()

            data = reader.GetOutput()
            bounds = data.GetBounds()
            xmin, xmax, ymin, ymax, zmin, zmax = bounds
            self.size = [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        except FileNotFoundError:
            print(f"Error: Volume file '{self.volume_file}' not found.")
            sys.exit(1)
        except ValueError:
            print("Error: volume.vtk file format is incorrect.")
            sys.exit(1)

    def read_step(self):
        if self.nsteps == "0":
            file_name = "elem.vtk"
        else:
            file_name = f"elem_{self.step:04}.vtk"

        try:
            with open(file_name, "r") as f:
                data_lines = f.readlines()

            self.t = float(data_lines[1].split()[0])

            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName(file_name)
            reader.Update()

            data = reader.GetOutput()

            self.npoints = data.GetNumberOfPoints()
            self.points = [data.GetPoint(i) for i in range(self.npoints)]
            self.ncells = data.GetNumberOfCells()
            self.cells = []
            for i in range(self.ncells):
                cell = data.GetCell(i)
                point_ids = cell.GetPointIds()
                self.cells.append([point_ids.GetId(0), point_ids.GetId(1)])
            point_data = data.GetPointData()
            scalars = point_data.GetScalars()
            self.point_scalars = [
                scalars.GetValue(i) / 2.0 for i in range(self.npoints)
            ]
        except FileNotFoundError:
            print(f"Error: Step file '{file_name}' not found.")
            sys.exit(1)
        except ValueError:
            print(f"Error: '{file_name}' file format is incorrect.")
            sys.exit(1)

    def set_up_plot(self):
        ax.clear()

        dx, dy, dz = (
            [self.size[0][0], self.size[1][0]],
            [self.size[0][1], self.size[1][1]],
            [self.size[0][2], self.size[1][2]],
        )
        aspect_ratios = [
            (self.size[1][i] - self.size[0][i]) / (self.size[1][0] - self.size[0][0])
            for i in range(3)
        ]

        ax.set_box_aspect(aspect_ratios)
        ax.set_xlim([dx[0], dx[1]])
        ax.set_ylim([dy[0], dy[1]])
        ax.set_zlim([dz[0], dz[1]])
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
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

    def do_plot_lines(self):
        for cell in self.cells:
            p0, p1 = self.points[cell[0]], self.points[cell[1]]

            if not self.plot_outside:
                mid_point = [(p0[i] + p1[i]) * 0.5 for i in range(3)]
                if any(
                    mid_point[i] < self.size[0][i] or mid_point[i] > self.size[1][i]
                    for i in range(3)
                ):
                    continue
                if any(
                    abs(p1[i] - p0[i]) > 0.5 * abs(self.size[1][i] - self.size[0][i])
                    for i in range(3)
                ):
                    continue
            dx, dy, dz = [p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]]
            ax.plot(dx, dy, dz, color=self.line_color, lw=self.line_width)

    def do_plot_points(self):
        px, py, pz = zip(
            *[
                p
                for p in self.points
                if self.plot_outside
                or all(self.size[0][i] <= p[i] <= self.size[1][i] for i in range(3))
            ]
        )
        ps = [
            self.point_scalars[i]
            for i in range(self.npoints)
            if self.plot_outside
            or all(
                self.size[0][j] <= self.points[i][j] <= self.size[1][j]
                for j in range(3)
            )
        ]
        ax.scatter(
            px,
            py,
            pz,
            s=self.point_size,
            c=ps,
            cmap="rainbow",
            vmin=0.0,
            vmax=1.0,
            alpha=1,
            edgecolor=self.point_edge_color,
        )

    def do_plot_step_info(self):
        t_num = f"{self.t:.5e}".split("e")[0]
        t_exp = int(f"{self.t:.5e}".split("e")[1])
        t_str = f"Step: {self.step}/{self.nsteps}, Time: {t_num} x10$^{{{t_exp}}}$ (sec)\n(Nodes: {self.npoints}, Elements: {self.ncells})"
        ax.set_title(t_str)

    def plot_step(self):
        self.set_up_plot()
        if self.plot_points:
            self.do_plot_points()
        if self.plot_lines:
            self.do_plot_lines()
        self.do_plot_step_info()
        plt.draw()

    def make_images_folder(self):
        if not os.path.isdir("images"):
            os.mkdir("images")

    def save_image(self):
        file_name = f"images/image_{self.step:04}.jpeg"
        plt.savefig(file_name)
        print(f"Step {self.step}/{self.nsteps} was saved: {file_name}")

    def save_images(self):
        for i in range(self.nsteps):
            self.step = i + 1
            self.read_step()
            self.plot_step()
            self.save_image()

    def make_movie(self):
        file_list = glob.glob(r"images/image_*.jpeg")
        file_list.sort()

        img = cv2.imread(file_list[0])
        h, w, channels = img.shape[:3]

        codec = cv2.VideoWriter_fourcc(*'mp4v')
        writer = cv2.VideoWriter("movie.mp4", codec, self.fps, (w, h), 1)
        bar = tqdm(total=len(file_list), dynamic_ncols=True)
        for f in file_list:
            img = cv2.imread(f)
            writer.write(img)
            bar.update(1)
        bar.close()
        writer.release()

    def on_key(self, event):
        if event.key == "up":
            self.step = min(self.step + 1, self.nsteps)
            self.read_step()
        elif event.key == "down":
            self.step = max(self.step - 1, 1)
            self.read_step()
        elif event.key == "n":
            self.plot_points = not self.plot_points
        elif event.key == "l":
            self.plot_lines = not self.plot_lines
        elif event.key == "a":
            self.plot_outside = not self.plot_outside
        elif event.key == "i":
            self.make_images_folder()
            self.save_image()
            return
        elif event.key == "m":
            self.make_images_folder()
            self.save_images()
            self.make_movie()
            self.step = 1
            self.read_step()
        elif event.key == "q":
            plt.close(fig)
            print("")
            return
        else:
            return
        self.plot_step()


def number_of_steps() -> int:
    nsteps = 1
    while True:
        file_name = f"elem_{nsteps:04}.vtk"
        if not os.path.isfile(file_name):
            break
        nsteps += 1
    if nsteps == 1 and not os.path.isfile("elem.vtk"):
        print("ERROR: Cannot find any input file.")
        sys.exit(1)
    return nsteps - 1


def main():
    element_plot = ElementPlotter(number_of_steps())
    element_plot.read_volume()
    element_plot.read_step()
    element_plot.plot_step()

    print("")
    print("Usage:")
    print("- Up key: plot the next step")
    print("- Down key: Plot the previous step")
    print("- n key: Plot/Not plot nodal points")
    print("- l key: Plot/Not plot element lines")
    print("- a key: Plot/Not plot points and lines outside the volume")
    print("- i key: Save figure as image file")
    print("- m key: Make movie file")
    print("- q key: Exit")
    print("")

    fig.canvas.mpl_connect("key_press_event", element_plot.on_key)
    plt.show()


if __name__ == "__main__":
    main()
