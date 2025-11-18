import numpy as np #type: ignore
import matplotlib.pyplot as plt#type: ignore
import sys

cutoff_radios_list = np.array(list(range(3,31)))
force_dipole = []

for n in cutoff_radios_list:
    filename = f"stress{n}"
    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps_related/{filename}.txt",'r') as f:
        lines = f.readlines()
        tokens = lines[1].split()
        force_dipole.append(float(tokens[0])) 

cutoff_radios_list = cutoff_radios_list*1.0e-10

plt.scatter(cutoff_radios_list, force_dipole, marker=".", s = 20)
plt.xlabel("Distance (x-x')[m]")
plt.ylabel("force dipole")
plt.show()