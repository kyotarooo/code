# =============================================================
# mbplot.py: Dislocation density plotter for Element-based PDD
# -------------------------------------------------------------
# Usage:
# - "l": Plot/Not plot lines
# - "n": Plot/Not plot points
# - "i": Save figure as image file
# - "q": Exit
# -------------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# =============================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Please modify as you like
LINE_WIDTH = 1.5
LINE_COLOR = "black"
POINT_SIZE = 12
POINT_COLOR = "black"
POINT_EDGE_COLOR = "black"
FILE_NAME = "density.out"


fig = plt.figure("Element-based PDD - Dislocation Density -")
ax = fig.add_subplot()


class DislocationDensity:
    def __init__(
        self,
        line_width=LINE_WIDTH,
        line_color=LINE_COLOR,
        point_size=POINT_SIZE,
        point_color=POINT_COLOR,
        point_edge_color=POINT_EDGE_COLOR,
        file_name=FILE_NAME,
    ):
        self.line_width = line_width
        self.line_color = line_color
        self.point_size = point_size
        self.point_color = point_color
        self.point_edge_color = point_edge_color
        self.file_name = file_name
        self.plot_lines = True
        self.plot_points = True

        self.time_data = []
        self.density_data = []

    def read_data(self):
        try:
            with open(self.file_name, "r") as f:
                for line in f:
                    t, density = map(float, line.split())
                    self.time_data.append(t)
                    self.density_data.append(density)
        except FileNotFoundError:
            print(f"Error: Data file '{self.file_name}' not found.")
        except ValueError as e:
            print(f"Error reading data: {e}")

    def plot(self):
        ax.clear()

        if self.plot_lines:
            ax.plot(
                self.time_data,
                self.density_data,
                color=self.line_color,
                lw=self.line_width,
            )
        if self.plot_points:
            ax.scatter(
                self.time_data,
                self.density_data,
                color=self.point_color,
                s=self.point_size,
                edgecolor=self.point_edge_color,
            )

        ax.grid(linestyle="--", linewidth=0.5)
        ax.set_xlabel("Time (sec)")
        ax.set_ylabel("Dislocation density (m$^{-2}$)")
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0)
        plt.draw()

    def save_image(self):
        plt.savefig("density.png")
        
    def on_key(self, event):
        if event.key == "l":
            self.plot_lines = not self.plot_lines
        elif event.key == "n":
            self.plot_points = not self.plot_points
        elif event.key == "i":
            self.save_image()
        elif event.key == "q":
            plt.close(fig)
        else:
            return
        self.plot()


def main():
    dd = DislocationDensity()
    dd.read_data()
    print("Usage:")
    print("- l: Plot/Not plot lines")
    print("- n: Plot/Not plot points")
    print("- i: Save figure as image file")
    print("- q: Exit")
    dd.plot()
    fig.canvas.mpl_connect("key_press_event", dd.on_key)
    plt.show()


if __name__ == "__main__":
    main()
