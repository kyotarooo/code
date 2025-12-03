# ======================================================================================
# loop_info.py
# --------------------------------------------------------------------------------------
# This python script extracts the dislocation loop information from elem_****.vtk files.
# The dislocation loop information is:
# 1. Number of open loops
# 2. Number of close loops
# 3. Length of each dislocation loop
# 4. Loop ID for each dislocation element
# The results are stored in loop_****.txt files.
# ---------------------------------------------------------------------------------------
# Developer: A. Takahashi (Tokyo University of Science)
# =======================================================================================

import vtk
import os
import math
from collections import defaultdict


class LoopData:
    """
    Stores element-node-geometry information for loop reconstruction.
    """

    def __init__(self):
        self.coords = []  # list of (x,y,z)
        self.nodes = []  # element -> [node0, node1]
        self.elements = []  # element -> [prev_element, next_element]
        self.loops = []  # element -> loop_id
        self.types = []  # loop_id -> "open" / "close"

        self.n_nodes = 0
        self.n_elements = 0
        self.n_loops = 0
        self.n_open_loops = 0
        self.n_close_loops = 0

        # NEW: reverse map from node → elements attached
        self.node_to_elements = defaultdict(list)

    def build_node_map(self):
        """Build fast node → element adjacency dictionary."""
        for eid, (n0, n1) in enumerate(self.nodes):
            self.node_to_elements[n0].append((eid, 0))  # element, position 0
            self.node_to_elements[n1].append((eid, 1))  # element, position 1


def find_start_element(loop_data):
    """
    Find an element whose loop ID is unset (-1), and is a good starting node.
    """
    for eid in range(loop_data.n_elements):
        if loop_data.loops[eid] == -1:
            nxt = loop_data.elements[eid][0]
            visited = set()

            while nxt != -1 and nxt not in visited:
                visited.add(nxt)
                nxt = loop_data.elements[nxt][0]

            return eid

    return -1


def make_loops(loop_data):
    loop_data.loops = [-1] * loop_data.n_elements
    loop_data.types = []

    loop_id = 0
    open_ct = 0
    close_ct = 0

    while True:
        start = find_start_element(loop_data)
        if start == -1:
            break

        eid = start
        while True:
            loop_data.loops[eid] = loop_id
            nxt = loop_data.elements[eid][1]

            if nxt == -1:
                loop_data.types.append("open")
                open_ct += 1
                break
            elif nxt == start:
                loop_data.types.append("close")
                close_ct += 1
                break

            eid = nxt

        loop_id += 1

    loop_data.n_loops = open_ct + close_ct
    loop_data.n_open_loops = open_ct
    loop_data.n_close_loops = close_ct


def loop_length(loop_data, loop_id):
    """Compute geometric length of a given loop."""
    length = 0.0

    for eid in range(loop_data.n_elements):
        if loop_data.loops[eid] == loop_id:
            n0, n1 = loop_data.nodes[eid]
            x0 = loop_data.coords[n0]
            x1 = loop_data.coords[n1]
            length += (x1[0] - x0[0]) ** 2 + (x1[1] - x0[1]) ** 2 + (x1[2] - x0[2]) ** 2

    return math.sqrt(length)


def number_of_steps():
    i = 1
    while True:
        if not os.path.isfile(f"elem_{i:04}.vtk"):
            return i - 1
        i += 1


def find_element_fast(loop_data, node_id, pos):
    """
    Return the element connected to (node_id) on side (pos).
    Uses precomputed dictionary for O(1) access.
    """
    for eid, p in loop_data.node_to_elements[node_id]:
        if p == pos:
            return eid
    return -1


def main():
    n_steps = number_of_steps()

    for step in range(1, n_steps + 1):

        # --- Read VTK data ---
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(f"elem_{step:04}.vtk")
        reader.Update()

        vtk_data = reader.GetOutput()

        loop_data = LoopData()
        loop_data.n_nodes = vtk_data.GetNumberOfPoints()
        loop_data.n_elements = vtk_data.GetNumberOfCells()

        # --- Read coordinates ---
        loop_data.coords = [vtk_data.GetPoint(i) for i in range(loop_data.n_nodes)]

        # --- Read elements (line segments) ---
        for eid in range(loop_data.n_elements):
            cell = vtk_data.GetCell(eid)
            pid = cell.GetPointIds()
            loop_data.nodes.append([pid.GetId(0), pid.GetId(1)])

        # Build node→element map for fast searching
        loop_data.build_node_map()

        # --- Build connectivity (prev, next) ---
        for eid in range(loop_data.n_elements):
            n0, n1 = loop_data.nodes[eid]

            prev_eid = find_element_fast(loop_data, n0, 1)
            next_eid = find_element_fast(loop_data, n1, 0)

            loop_data.elements.append([prev_eid, next_eid])
            loop_data.loops.append(-1)

        # --- Build loops ---
        make_loops(loop_data)

        # --- Output results ---
        with open(f"loop_{step:04}.txt", "w") as f:
            f.write(f"NOpenLoops: {loop_data.n_open_loops}\n")
            for lid in range(loop_data.n_loops):
                if loop_data.types[lid] == "open":
                    f.write(f"LoopID: {lid} Length: {loop_length(loop_data, lid):e}\n")

            f.write(f"NCloseLoops: {loop_data.n_close_loops}\n")
            for lid in range(loop_data.n_loops):
                if loop_data.types[lid] == "close":
                    f.write(f"LoopID: {lid} Length: {loop_length(loop_data, lid):e}\n")

            f.write(f"{loop_data.n_elements}\n")
            for eid in range(loop_data.n_elements):
                f.write(f"{loop_data.loops[eid]}\n")


if __name__ == "__main__":
    main()
