# =========================================================================================
# material.inp作成
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

############ Generate supercell points along diagonal ############
with open(f"{output_dir}/supercell.txt", "r") as f_cell:
    lines = f_cell.readlines()
    cell = np.array([float(v) for v in lines[5].split()])
    

def main():
    for cutoff in range(cutoff_num):
        mk_dir = f"/Users/kyou/Library/CloudStorage/Box-Box/code/InclusionStress-20251127/input2/disp_{cutoff}_{a_core_width}"
        os.makedirs(mk_dir, exist_ok=True)
        with open(f'{mk_dir}/material.inp', "w") as f_incl:
            f_incl.write("1.0 0.0 0.0\n")
            f_incl.write("0.0 1.0 0.0\n")
            f_incl.write("0.0 0.0 1.0\n")
            f_incl.write("1 1 1\n")
            f_incl.write(f"{200 * conv_ang_to_m} {200 * conv_ang_to_m} {200 * conv_ang_to_m}\n")
            f_incl.write("1.92965e+11 2.01281e-01")
            
    mk_dir = f"/Users/kyou/Library/CloudStorage/Box-Box/code/InclusionStress-20251127/input2/residual_{a_core_width}"
    os.makedirs(mk_dir, exist_ok=True)
    with open(f'{mk_dir}/material.inp', "w") as f_re:
        f_re.write("1.0 0.0 0.0\n")
        f_re.write("0.0 1.0 0.0\n")
        f_re.write("0.0 0.0 1.0\n")
        f_re.write("1 1 1\n")
        f_re.write(f"{200 * conv_ang_to_m} {200 * conv_ang_to_m} {200 * conv_ang_to_m}\n")
        f_re.write("1.92965e+11 2.01281e-01")
            
    for cutoff in range(cutoff_num):
        mk_dir = f"/Users/kyou/Library/CloudStorage/Box-Box/code/InclusionStress-20251127/input2/disp_{cutoff}_{a_core_width}"
        os.makedirs(mk_dir, exist_ok=True)        
        with open(f'{mk_dir}/grid.inp', "w") as f_incl:
            f_incl.write("2\n")
            f_incl.write("100 100 100\n")
            f_incl.write("10 10 10\n")
            
    mk_dir = f"/Users/kyou/Library/CloudStorage/Box-Box/code/InclusionStress-20251127/input2/residual_{a_core_width}"
    os.makedirs(mk_dir, exist_ok=True)           
    with open(f'{mk_dir}/grid.inp', "w") as f_re:
        f_re.write("2\n")
        f_re.write("100 100 100\n")
        f_re.write("10 10 10\n")
    
if __name__ == "__main__":
     main()


        
             

             

             


             

            