import numpy as np#type: ignore
import sys

def main():
    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_J/perfect_J_sumstress.txt", "r") as f_pe_J,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_eV/perfect_eV_sumstress.txt", "r") as f_pe_eV,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_J/defect_J_sumstress.txt", "r") as f_de_J,\
         open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_eV/defect_eV_sumstress.txt", "r") as f_de_eV:
        lines_perfect_J = f_pe_J.readlines()
        lines_perfect_eV = f_pe_eV.readlines()
        lines_defect_J = f_de_J.readlines()
        lines_defect_eV = f_de_eV.readlines()
        
        n = len(lines_perfect_J)

        print (len(lines_perfect_J))
        print (len(lines_defect_J))
        print (len(lines_perfect_eV))
        print (len(lines_defect_eV))
        p_output_eV = np.zeros((100, 6))
        p_output_J = np.zeros((100, 6))

        id_max = 0
        for i in range(n):
            id = int(lines_perfect_J[i].split()[0])
            id_max = max(id_max, id)

            for j in range(6):
                p_J = float(lines_perfect_J[i].split()[j + 1])
                d_J = float(lines_defect_J[i].split()[j + 1])
                p_eV = float(lines_perfect_eV[i].split()[j + 1])
                d_eV = float(lines_defect_eV[i].split()[j + 1])

                p_output_eV[id][j] = -(d_eV - p_eV)
                p_output_J[id][j] = -(d_J - p_J)

    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_eV/force_dipole_eV.txt", "w") as f_eV:
        for i in range(id_max + 1):
            for k in range(6):
                f_eV.write(f"{p_output_eV[i][k]} ")
            f_eV.write("\n")
    
    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_J/force_dipole_J.txt", "w") as f_J:
        for i in range(id_max + 1):
            for k in range(6):
                f_J.write(f"{p_output_J[i][k]} ")
            f_J.write("\n")



if __name__ == "__main__":
    main()















