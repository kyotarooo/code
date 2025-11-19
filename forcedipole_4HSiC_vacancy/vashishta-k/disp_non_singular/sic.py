#last updated: Apr 26, 2022

#!/usr/bin/env python

import numpy as np # type: ignore
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")


import math
from numpy import array #type: ignore
import numpy as np #type: ignore

class SiC:
    d = 3.07659964
    a1 = d
    a2 = math.sqrt(3) * d
    a3 = math.sqrt(6) / 3.0 * d

    def __init__(self, pattern, nx, ny, nz=None, lay_sum=None):
        pattern_index = {"a": 0, "b": 1, "c": 2}
        self.pattern = array([pattern_index[x] for x in pattern])

        # Define box
        self.nx = nx
        self.ny = ny
        self.lay_num = len(pattern)
        self.lx = nx * SiC.a1
        self.ly = ny * SiC.a2
        self.type = ["Si", "C"]

        if lay_sum is None and nz is not None:
            self.nz = nz * self.lay_num
        elif lay_sum is not None:
            self.nz = lay_sum
        else:
            raise ValueError("Either 'nz' or 'lay_sum' is required")

        self.lz = self.nz * SiC.a3


        # Set up atoms
        self.atom_arr_list = []

        x = array([SiC.a1, 0, 0])
        y = array([0, SiC.a2, 0])
        z = array([0, 0, SiC.a3])
        o = array([0, 0, 0])

        lay_a = array([o, (x + y) / 2])
        lay_b = array([y / 3, x / 2 + 5 * y / 6])
        lay_c = array([2 * y / 3, x / 2 + y / 6])
        lays = array([lay_a, lay_b, lay_c])

        tmp = [self.pattern, np.hstack([self.pattern[1:], self.pattern[0:1]])]
        unit = [lays[i].reshape(-1, 3) for i in tmp]

        for n, u in enumerate(unit):
            if self.nz % len(pattern) == 0:
                arry = np.fromfunction(lambda i, j, k, n: x * i + y * j + z * (k // 2),
                                       (self.nx, self.ny, 2 * self.nz, 3))
                arry = arry.reshape(-1, *(u.shape))
                arry += u
            else:
                arry = np.fromfunction(lambda i, j, k, n: x * i + y * j + z * (k // 2),
                                       (self.nx, self.ny, 2 * self.nz, 3))

                tmp = np.ones((2 * (self.nz // len(pattern) + 1) * len(pattern), 3))
                tmp = tmp.reshape(-1, *(u.shape))
                u = tmp * u
                u = u.reshape(-1, 3)[:2 * self.nz]

                arry = arry.reshape(-1, *(u.shape))
                arry += u
                
            self.atom_arr_list.append(arry.reshape(-1, 3) + z * (n * 0.25))
            
        self.atomnum = 0
        for arry in self.atom_arr_list:
            self.atomnum += arry.shape[0]

    def setBox(self, Lx=None, Ly=None, Lz=None):
        if Lx != None:
            self.lx = Lx
        if Ly != None:
            self.ly = Ly
        if Lz != None:
            self.lz = Lz

    def rev(self):
        """Change Si and C atomlist"""
        
        buf = self.atom_arr_list[0]
        self.atom_arr_list[0] = self.atom_arr_list[1]
        self.atom_arr_list[1] = buf

    def outAtomeye(self, f):
        #print("Number of particles =", self.atomnum, file=f)
        print("Number of particles =", self.atom_arr_list[0].shape[0]+self.atom_arr_list[1].shape[0], file=f)
        print("A = 1 Angstrom (basic length-scale)", file=f)
        print("H0(1,1) =", self.lx, "A", file=f)
        print("H0(1,2) = 0 A", file=f)
        print("H0(1,3) = 0 A", file=f)
        print("H0(2,1) = 0 A", file=f)
        print("H0(2,2) =", self.ly, "A", file=f)
        print("H0(2,3) = 0 A", file=f)
        print("H0(3,1) = 0 A", file=f)
        print("H0(3,2) = 0 A", file=f)
        print("H0(3,3) =", self.lz, "A", file=f)
        print(".NO_VELOCITY.", file=f)
        print("entry_count = 4", file=f)
        print("auxiliary[0] = id", file=f)
        count = 1
        for i, arry in enumerate(self.atom_arr_list):
            for atom in arry:
                print("%d" % (i + 1), file=f)
                print(self.type[i], file=f)
                print("%f %f %f %d" % (atom[0] / self.lx, atom[1] / self.ly, atom[2] / self.lz, count), file=f)
                count += 1

    def get_atom_num(self):
        self.atomnum = 0
        for arry in self.atom_arr_list:
            self.atomnum += arry.shape[0]
        return self.atomnum

    def outCoords(self):
        print(self.get_atom_num())
        count = 1
        for i, arry in enumerate(self.atom_arr_list):
            for atom in arry:
                print("%d %f %f %f" % (count, atom[0], atom[1], atom[2]))
                count += 1

    def outLammps(self, fout):
        for i, arry in enumerate(self.atom_arr_list):
            for atom in arry:
                print("create_atoms %d single %f %f %f remap yes units box" % (
                    i + 1, atom[0], atom[1], atom[2]), file=fout)

    def outLammpsd(self, fout=None):
        print("LAMMPS data file via write_data, version 30 Jul 2016, timestep = 0", file=fout)
        print("%d atoms" % self.get_atom_num(), file=fout)
        print("2 atom types\n", file=fout)
        print("0.0 %f xlo xhi" % self.lx, file=fout)
        print("0.0 %f ylo yhi" % self.ly, file=fout)
        print("0.0 %f zlo zhi" % self.lz, file=fout)
        # print("0.0 0.0 0.0 xy xz yz", file=fout)
        print("\nMasses\n", file=fout)
        print("1 28.0855", file=fout)
        print("2 12.0107", file=fout)
        print("\nAtoms\n", file=fout)
        
        count = 1
        for i, arry in enumerate(self.atom_arr_list):
            for atom in arry:
                print("%d %d %f %f %f" % (count, i+1, atom[0], atom[1], atom[2]), file=fout)
                count += 1

    def outLammpsdq(self, fout=None):
        q = [0, 0]
        print("LAMMPS data file via write_data, version 30 Jul 2016, timestep = 0", file=fout)
        print("%d atoms" % self.get_atom_num(), file=fout)
        print("2 atom types\n", file=fout)
        print("0.0 %f xlo xhi" % self.lx, file=fout)
        print("0.0 %f ylo yhi" % self.ly, file=fout)
        print("0.0 %f zlo zhi" % self.lz, file=fout)
        # print("0.0 0.0 0.0 xy xz yz", file=fout)
        print("\nMasses\n", file=fout)
        print("1 28.0855", file=fout)
        print("2 12.0107", file=fout)
        print("\nAtoms\n", file=fout)
        
        count = 1
        for i, arry in enumerate(self.atom_arr_list):
            for atom in arry:
                print("%d %d %f %f %f %f" % (count, i+1, q[i], atom[0], atom[1], atom[2]), file=fout)
                count += 1


#n_list = [[1, 1, 1], [2, 1, 1], [2, 2, 1], [3, 2, 1], [2, 3, 2], [4, 3, 2], [3, 3, 3], [5, 3, 3], [3, 4, 3], [5, 4, 3], [4, 6, 3], [8, 6, 3], [5, 7, 3], [10, 7, 3]]
n_list = [[1, 1, 1], [2, 1, 1], [3, 2, 1], [4, 3, 2], [5, 3, 3], [5, 4, 3], [8, 6, 3], [10, 7, 3]]
if __name__ == "__main__":
    id = 0
    with open(f"{output_dir}/4hsic_q/filename.txt" , "w") as f_f:
        for nx, ny, nz in n_list:
            # ------------config--------------
            pattern = "abac"
            # --------------------------------
            a = SiC(pattern, nx, ny, nz)
            id = id + 1
            a.outAtomeye(open("test0.cfg", "w"))
            a.outLammpsd(open(f"{output_dir}/4hsic/{id}_in.sic", "w"))
            a.outLammpsdq(open(f"{output_dir}/4hsic_q/{id}_in.sicq", "w"))
            f_f.write(f"{id}_in.sicq\n")
            
            with open(f"{output_dir}/include/{id}_include_n", "w") as f:
                print("variable nx equal %d" % nx, file=f)
                print("variable ny equal %d" % ny, file=f)
                print("variable nz equal %d" % nz, file=f)
