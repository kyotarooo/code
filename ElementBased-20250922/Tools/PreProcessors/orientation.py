import numpy as np
import TkEasyGUI as eg

layout = [
    [
        eg.Text("X axis"),
        eg.Text("X"),
        eg.Input("1.0", key="-XAXISX-"),
        eg.Text("Y"),
        eg.Input("0.0", key="-XAXISY-"),
        eg.Text("Z"),
        eg.Input("0.0", key="-XAXISZ-"),
    ],
    [
        eg.Text("Y axis"),
        eg.Text("X"),
        eg.Input("0.0", key="-YAXISX-"),
        eg.Text("Y"),
        eg.Input("1.0", key="-YAXISY-"),
        eg.Text("Z"),
        eg.Input("0.0", key="-YAXISZ-"),
    ],
    [
        eg.Text("Z axis"),
        eg.Text("X"),
        eg.Input("0.0", key="-ZAXISX-"),
        eg.Text("Y"),
        eg.Input("0.0", key="-ZAXISY-"),
        eg.Text("Z"),
        eg.Input("1.0", key="-ZAXISZ-"),
    ],
    [
        eg.Text("Target vector"),
        eg.Text("X"),
        eg.Input("0.0", key="-X-"),
        eg.Text("Y"),
        eg.Input("0.0", key="-Y-"),
        eg.Text("Z"),
        eg.Input("0.0", key="-Z-"),
    ],
    [eg.HSeparator()],
    [eg.Button("OK"), eg.Button("Cancel")],
]

def calculate_vector(values):
    v = np.array([float(values["-XAXISX-"]), float(values["-XAXISY-"]), float(values["-XAXISZ-"])])
    ax = v / np.linalg.norm(v)
    v = np.array([float(values["-YAXISX-"]), float(values["-YAXISY-"]), float(values["-YAXISZ-"])])
    ay = v / np.linalg.norm(v)
    v = np.array([float(values["-ZAXISX-"]), float(values["-ZAXISY-"]), float(values["-ZAXISZ-"])])
    az = v / np.linalg.norm(v)
    T = [ax, ay, az]

    x = [float(values["-X-"]), float(values["-Y-"]), float(values["-Z-"])]
    xp = []
    for i in range(3):
        v = 0.0
        for j in range(3):
            v += T[i][j] * x[j]
        xp.append(v)
    v = np.array(xp)
    e = v / np.linalg.norm(v)
    print(f"{xp[0]:e} {xp[1]:e} {xp[2]:e}")
    print(f"({e[0]:e} {e[1]:e} {e[2]:e})")


def main():
    window = eg.Window("Crystal Orientation Calculator", layout=layout)

    while window.is_alive():
        event, values = window.read()
        if event == "OK":
            calculate_vector(values)
        elif event == "Cancel":
            break
    window.close()

if __name__ == "__main__":
    main()
