import numpy as np


class ElemData:
    def __init__(self):
        self.n_nodes = 0
        self.types = []
        self.coords = []
        self.n_elements = 0
        self.nodes = []
        self.burgers_vecs = []
        self.slip_vecs = []

    def Read(self):
        with open("elem.inp", "r") as f:
            lines = f.readlines()

        index = 0
        self.n_nodes = int(lines[index].split()[0])
        index += 1
        for i in range(self.n_nodes):
            self.types.append(int(lines[index].split()[0]))
            self.coords.append([float(lines[index].split()[j + 1]) for j in range(3)])
            index += 1

        self.n_elements = int(lines[index].split()[0])
        index += 1
        for i in range(self.n_elements):
            self.nodes.append([int(lines[index].split()[j]) for j in range(2)])
            self.burgers_vecs.append(
                [float(lines[index].split()[j + 2]) for j in range(3)]
            )
            self.slip_vecs.append(
                [float(lines[index].split()[j + 5]) for j in range(3)]
            )
            index += 1

    def CheckElementOnSlipPlane(self):
        print("Checking elements are on slip planes.")
        print(
            "If the output (dot-product of tangent vector and slip plane normal vector) is zero,"
        )
        print("the element is properly placed on the slip plane.")
        for i in range(self.n_elements):
            node_id0 = self.nodes[i][0]
            node_id1 = self.nodes[i][1]
            coord0 = np.array(self.coords[node_id0])
            coord1 = np.array(self.coords[node_id1])
            s = np.array(self.slip_vecs[i])

            t = coord1 - coord0
            dot_product = np.dot(t, s)

            print(f"Element({i}): {dot_product}")

    def CheckBurgersVectorOnSlipPlane(self):
        print("Checking Burgers vectors are on slip planes.")
        print(
            "If the output (dot-product of Burgers vector and slip plane normal vector) is zero,"
        )
        print("the Burgers vector is on the slip plane.")
        for i in range(self.n_elements):
            b = np.array(self.burgers_vecs[i])
            s = np.array(self.slip_vecs[i])

            dot_product = np.dot(b, s)

            print(f"Element({i}): {dot_product}")


def main():
    elem_data = ElemData()

    elem_data.Read()
    print("")
    print("[Step 1]")
    elem_data.CheckElementOnSlipPlane()
    print("")
    print("[Step 2]")
    elem_data.CheckBurgersVectorOnSlipPlane()


if __name__ == "__main__":
    main()
