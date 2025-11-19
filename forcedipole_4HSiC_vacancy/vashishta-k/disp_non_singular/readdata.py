import numpy as np#type: ignore
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

######## 削除する原子の指定 ########  どの位置の原子を削除するか ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
atom_type_to_delete = int(os.environ.get("ATOM"))
if atom_type_to_delete is None:
    raise ValueError("環境変数 ATOM error")

# ####### 定数の定義(等方性材料) ############
# E = 70.6e+09  # pa
# nu = 0.3 
# g = E / (2 * (1 + nu))  

######### 定数の定義(異方性材料) ############
c_11 = 522 # Gpa
c_12 = 127 # Gpa
c_44 = 190 # Gpa

c_11 = c_11 * 1.0e+9 #pa
c_12 = c_12 * 1.0e+9 #pa
c_44 = c_44 * 1.0e+9 #pa


#======== 単位 ========
# 1 bar = 10^5 pa

######## Voigt average ########
lambda_v = (c_11 + 4 * c_12 - 2 * c_44) / 5
mu_v = (c_11 - c_12 + 3 * c_44) / 5

######## Reuss average ########
lambda_r = (c_11 * c_11 + c_11 * c_12 - 2 * c_11 * c_44 - 2 * c_12 * c_12 + 6 * c_12 * c_44) / (3 * (c_11 - c_12) + 4 * c_44)
mu_r = (5 * c_44 * (c_11 - c_12)) / (3 * (c_11 - c_12) + 4 * c_44)

######## ヤング率, ポアソン比, 剛性率　の算出 ########
E_v = mu_v * (3 * lambda_v + 2 * mu_v) / (lambda_v + mu_v)
E_r = mu_r * (3 * lambda_r + 2 * mu_r) / (lambda_r + mu_r)
nu_v = lambda_v / (2 * (lambda_v + mu_v))
nu_r = lambda_r / (2 * (lambda_r + mu_r))
g_v = E_v / (2 * (1 + nu_v))
g_r = E_r / (2 * (1 + nu_r))

E = (E_v + E_r) / 2
nu = (nu_v + nu_r) / 2
g = (g_v + g_r) / 2

print(f"{'lambda:':<10}{'Voigt average':<20}={E_v / 1e9:10.3f} GPa")
print(f"{'lambda:':<10}{'Reuss average':<20}={E_r / 1e9:10.3f} GPa")
print(f"{'nu:':<10}{'Voigt average':<20}={nu_v:10.3f} [-]")
print(f"{'nu:':<10}{'Reuss average':<20}={nu_r:10.3f} [-]")
print("******** input data ********")
print(f"g      =    {g:.5e} [pa]")
print(f"nu     =    {nu:.5e} [-]")

######## 単位変換 ########
conv_ang_to_m = 1.0E-10

######## dump(perfect)ファイルの読み込み関数 ########
def read_perfect_data(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    fast_num = 0
    for i,line in enumerate(lines):
        if "ITEM: ATOMS" in line:
            fast_num = i
            break

    atoms = {} 
    for line in lines[fast_num+1:]:
        if "ITEM" in line:
            break
        tokens = line.strip().split()
        atom_id = int(tokens[0])
        x,y,z = map(float, tokens[2:5])
        atoms[atom_id] = np.array([x,y,z])
    
    return atoms

######## dump(defect)ファイルの読み込み関数 ########
def read_defect_data(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    last_num = 0
    for i,line in enumerate(lines):
        if "ITEM: ATOMS" in line:
            last_num = i
        
    atoms = {}
    for line in lines[last_num+1:]:
        list = line.strip().split()
        atom_id = int(list[0])
        x,y,z = map(float, list[2:5])
        atoms[atom_id] = np.array([x,y,z])

    return atoms

######## dump(deleted atom)ファイルの読み込み関数 ########
def read_vacancydata(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        items = lines[0].strip().split()
        
    x,y,z = map(float, items[3:6])
    atoms = np.array([x,y,z])
    return atoms

######## ファイル数の読み込み ########
with open(f"{output_dir}/4hsic_q/filename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file_num = len(lines)

######## supercell sizeの読み込み ########
disp_threshold = []
supercell_size = []
with open(f"{output_dir}/supercell.txt", 'r') as f_size:
    lines = f_size.readlines()
    data_amount = len(lines)
    for i in range(data_amount-file_num, data_amount):
        size = np.array(lines[i].split(),dtype =float)
        supercell_size.append(size)
        half_size = size/2
        disp_threshold.append(half_size)

    
######## 座標データの読み込み ########
coordinate_of_defect = []
coordinate_of_perfect = []
coordinate_of_vacancy = []
atom_num = []
common_ids = []
for i in range(file_num):
    i += 1
    coordinate_p = read_perfect_data(f"{output_dir}/dump/perfect_lattice/perfect_lattice_{atom_type_to_delete}/perfect_lattice{i}.dump")
    coordinate_of_perfect.append(coordinate_p)
    coordinate_d = read_defect_data(f"{output_dir}/dump/defect_lattice/defect_lattice_{atom_type_to_delete}/defect_lattice{i}.dump")
    coordinate_of_defect.append(coordinate_d)
    coordinate_v = read_vacancydata(f"{output_dir}/dump/deleted_atom/deleted_atom_{atom_type_to_delete}/deleted_atom{i}")
    coordinate_of_vacancy.append(coordinate_v)
    atom_num.append(len(coordinate_d))

    common_id = sorted(set(coordinate_p))
    common_ids.append(common_id)  

######## displacementの出力 ########
for i in range(file_num):
    with open(f"{output_dir}/dump/displacement/displacement_{atom_type_to_delete}/displacement{i}.inp", "w") as f_disp,\
        open(f"{output_dir}/dump/displacement_ovit/displacement_ovit_{atom_type_to_delete}/displacement_ovit{i}.xyz", "w") as f_disp_ovit:
        
        # 拡張xyzのヘッダ行 <pecies: 元素記号, pos: 座標, disp: 変位ベクトル（これがOVITOでベクトル表示される）>
        f_disp_ovit.write(f"{atom_num[i]}\n")
        f_disp_ovit.write("Properties=species:S:1:pos:R:3:disp:R:3\n")

        f_disp.write(f"{g} {nu}\n")
        xc1,xc2,xc3 = coordinate_of_vacancy[i]
        f_disp.write(f"{xc1*conv_ang_to_m} {xc2*conv_ang_to_m} {xc3*conv_ang_to_m}\n")
        f_disp.write(f"{atom_num[i]}\n")

        for atom_id in common_ids[i]:
            x = coordinate_of_perfect[i][atom_id]
            u = coordinate_of_defect[i][atom_id] - coordinate_of_perfect[i][atom_id]
            for k in range(3):
                if abs(u[k]) > disp_threshold[i][k]:
                    if u[k] > 0:
                        u[k] -= supercell_size[i][k]
                    elif u[k] < 0:
                         u[k] += supercell_size[i][k]

            f_disp_ovit.write(f"Al {x[0]:.5f} {x[1]:.5f} {x[2]:.5f} {u[0]:.5f} {u[1]:.5f} {u[2]:.5f}\n")
            f_disp.write(f"{x[0]*conv_ang_to_m:.15e} {x[1]*conv_ang_to_m:.15e} {x[2]*conv_ang_to_m:.15e} {u[0]*conv_ang_to_m:.15e} {u[1]*conv_ang_to_m:.15e} {u[2]*conv_ang_to_m:.15e}\n")

            

