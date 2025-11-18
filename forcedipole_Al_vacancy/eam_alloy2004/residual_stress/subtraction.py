import numpy as np#type: ignore
import sys
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

def main():
    with open(f"{output_dir}/sumstress/perfect_J_sumstress.txt", "r") as f_pe_J,\
         open(f"{output_dir}/sumstress/perfect_eV_sumstress.txt", "r") as f_pe_eV,\
         open(f"{output_dir}/sumstress/defect_J_sumstress.txt", "r") as f_de_J,\
         open(f"{output_dir}/sumstress/defect_eV_sumstress.txt", "r") as f_de_eV:
        lines_perfect_J = f_pe_J.readlines()
        lines_perfect_eV = f_pe_eV.readlines()
        lines_defect_J = f_de_J.readlines()
        lines_defect_eV = f_de_eV.readlines()
        
        n = len(lines_defect_J)

        print (len(lines_perfect_J))
        print (len(lines_defect_J))
        print (len(lines_perfect_eV))
        print (len(lines_defect_eV))
        p_output_eV = np.zeros((100, 9))
        p_output_J = np.zeros((100, 9))

        id_max = 0
        for i in range(n):
            id = int(lines_perfect_J[i].split()[0])
            id_max = max(id_max, id)

            for j in range(9):
                p_J = float(lines_perfect_J[i].split()[j + 1])
                d_J = float(lines_defect_J[i].split()[j + 1])
                p_eV = float(lines_perfect_eV[i].split()[j + 1])
                d_eV = float(lines_defect_eV[i].split()[j + 1])

                p_output_eV[id][j] = -(d_eV - p_eV)
                p_output_J[id][j] = -(d_J - p_J)

    print(f"idmax = {id_max}")

    with open(f"{output_dir}/force_dipole/force_dipole_eV.txt", "w") as f_eV:
        for i in range(id_max):
            for k in range(9):
                f_eV.write(f"{p_output_eV[i][k]} ")
            f_eV.write("\n")
    
    with open(f"{output_dir}/force_dipole/force_dipole_J.txt", "w") as f_J:
        for i in range(id_max):
            for k in range(9):
                f_J.write(f"{p_output_J[i][k]} ")
            f_J.write("\n")



if __name__ == "__main__":
    main()















