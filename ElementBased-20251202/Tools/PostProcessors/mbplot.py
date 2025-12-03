# =============================================================
# mbplot.py: Mechanical behavior plotter for Element-based PDD
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

# Please modify as you like
LINE_WIDTH = 1.5
LINE_COLOR = "black"
POINT_SIZE = 12
POINT_COLOR = "red"
POINT_EDGE_COLOR = "black"
FILE_NAME = "mechanical_behavior.out"


fig = plt.figure("Element-based PDD - Mechanical Behavior -")
ax = fig.add_subplot()


class MechanicalBehavior:
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
        self.plotlines = True
        self.plotpoints = True

        self.strain = []
        self.plastic_strain = []
        self.stress = []

    def read_data(self):
        try:
            with open(self.file_name, "r") as f:
                for line in f:
                    t, strain, plastic_strain, stress = map(float, line.split())
                    self.strain.append(strain * 100.0)
                    self.plastic_strain.append(plastic_strain * 100.0)
                    self.stress.append(stress * 1.0e-06)
        except FileNotFoundError:
            print(f"Error: File '{self.file_name}' not found.")

    def plot(self):
        ax.clear()

        if self.plotlines:
            ax.plot(self.strain, self.stress, color=self.line_color, lw=self.line_width, zorder=1,)
        if self.plotpoints:
            ax.scatter(
                self.strain,
                self.stress,
                color=self.point_color,
                s=self.point_size,
                edgecolor=self.point_edge_color,
                zorder=2,
            )

        ax.grid(linestyle="--", linewidth=0.5)
        ax.set_xlabel("Strain (%)")
        ax.set_ylabel("Stress (MPa)")
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0)
        plt.draw()

    def save_image(self):
        plt.savefig("mechanical_behavior.png")
        
    def on_key(self, event):
        if event.key == "l":
            self.plotlines = not self.plotlines
        elif event.key == "n":
            self.plotpoints = not self.plotpoints
        elif event.key == "up":
            self.point_size += 1
        elif event.key == "down":
            self.point_size -= 1
        elif event.key == "i":
            self.save_image()
        elif event.key == "q":
            plt.close(fig)
            return
        else:
            return
        self.plot()


def main():
    mb = MechanicalBehavior()
    mb.read_data()
    print("Usage:")
    print("- l: Plot/Not plot lines")
    print("- n: Plot/Not plot points")
    print("- up: Increase the point size")
    print("- down: Decrease the point size")
    print("- i: Save figure as image file")
    print("- q: Exit")
    mb.plot()
    fig.canvas.mpl_connect("key_press_event", mb.on_key)
    plt.show()


if __name__ == "__main__":
    main()
