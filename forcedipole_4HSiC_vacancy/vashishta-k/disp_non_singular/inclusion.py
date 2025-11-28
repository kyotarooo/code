# =========================================================================================
# inclusion.inp作成
# =========================================================================================

import numpy as np #type: ignore
import matplotlib.pyplot as plt #type: ignore 
import os

############ Get working directory path ############
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

############ Get force dipole(residual stress) directory path ############
dipole_r_path = os.environ.get("DIPOLE_R_PATH")
if dipole_r_path is None:
    raise ValueError("環境変数 DIPOLE_R_PATH error")

############ Define material constant ################
lattice_const = os.environ.get("LATTICE_CONST") # ang
lattice_const = float(lattice_const) * 1.0e-10 # m


############# Unit conversion constants ############
conv_eV_to_J = 1.602176634E-19
conv_bars_to_Pa = 1.0E+5
conv_Vang_to_Vm = 1.0E-30
conv_barVang_to_NVm = conv_bars_to_Pa * conv_Vang_to_Vm
conv_J_to_eV = 1/conv_eV_to_J
conv_ang_to_m = 1.0E-10

######## 削除する原子の指定 ########  
atom_type_to_delete = os.environ.get("ATOM")
if atom_type_to_delete is None:
    raise ValueError("環境変数 ATOM error")

######## コア幅 ########  ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
a_core_width = float(os.environ.get("Burgers"))
if a_core_width is None:
    raise ValueError("環境変数 core error")


######## dump(deleted atom)ファイルの読み込み関数 ########
def read_vacancydata(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        items = lines[0].strip().split()
        
    x,y,z = map(float, items[3:6])
    atoms = np.array([x,y,z])
    return atoms
    
######### ファイル数の読み込み ########
with open(f"{output_dir}/4hsic_vacancy/atom_type_delete_{atom_type_to_delete}/filename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file = lines[0].split()
    file_num = len(file)

############ Load Force Dipole tensor (displacement) [J] ############
with open(f"{output_dir}/force_dipole/force_dipole_{atom_type_to_delete}/a_{a_core_width}/data_p_J.inp", "r") as f_P,\
    open(f"{output_dir}/force_dipole/mode.txt", "w") as f_mode:
    lines = f_P.readlines()
    cutoff_num = 4
    P = {i: [] for i in range(cutoff_num)} 

    for i in range(len(lines)):
        values = np.array([float(v) for v in lines[i].split()])
        mod = i % cutoff_num
        P[mod].append(values)

############ Load Force Dipole tensor (residual_stress) [J] ############
with open(f"{dipole_r_path}/force_dipole/force_dipole_{atom_type_to_delete}/force_dipole_J.txt", "r") as f_r:
    lines = f_r.readlines()
    P_r = []
    for i in range(len(lines)):
        values = np.array([float(v) for v in lines[i].split()])
        P_r.append(values)
    num_P_r = (len(P_r))

############ Load The number of atoms ############
with open(f"{output_dir}/N_atom_perfect.txt", "r") as f_atoms:
    atom = f_atoms.readlines()
    for i in range(len(atom)):
        atoms = [float(line.strip()) for line in atom]


def main():

    coordinate_v = read_vacancydata(f"{output_dir}/dump/deleted_atom/deleted_atom_{atom_type_to_delete}/deleted_atom{6}")  # ←←←←←←←←←←←←←←←←←←N = 959 deleted_atomは1~8
    xc1, xc2, xc3 = coordinate_v
    
    for cutoff in range(cutoff_num):
        mk_dir = f"/Users/kyou/Library/CloudStorage/Box-Box/code/InclusionStress-20251124/disp_{cutoff}_{a_core_width}"
        os.makedirs(mk_dir, exist_ok=True)
        with open(f'{mk_dir}/inclusion.inp', "w") as f_incl:
            f_incl.write("1\n")
            f_incl.write("elastic_dipole\n")
            f_incl.write("0.0 ")
            f_incl.write(f"{a_core_width:e}\n")
            f_incl.write(f"{xc1 * conv_ang_to_m} {xc2 * conv_ang_to_m} {xc3 * conv_ang_to_m}\n")
            
            for i in range(9):
                f_incl.write(f"{P[cutoff][5][i]} ")   # ←←←←←←←←←←←←←←←←←←N = 959
                if i == 2 or i == 5:
                    f_incl.write("\n")
            f_incl.write("\n")
            
    mk_dir = f"/Users/kyou/Library/CloudStorage/Box-Box/code/InclusionStress-20251124/residual_{a_core_width}"
    os.makedirs(mk_dir, exist_ok=True) 
    with open(f'{mk_dir}/inclusion.inp', "w") as f_re:
        f_re.write("1\n")
        f_re.write("elastic_dipole\n")
        f_re.write("0.0 ")
        f_re.write(f"{a_core_width:e}\n")
        f_re.write(f"{xc1 * conv_ang_to_m} {xc2 * conv_ang_to_m} {xc3 * conv_ang_to_m}\n")
        for i in range(9):
            f_re.write(f"{P_r[5][i]} ")
            if i == 2 or i == 5:
                    f_re.write("\n")

if __name__ == "__main__":
     main()


        
             

             

             


             

            