import numpy as np#type: ignore

######### 定数の定義 ############
E = 70.6e+09  # Pa ←←←材料力学　裳華房！！
nu = 0.3 
#g = E / (2 * (1 + nu))  # 剛性率(等方性)
g = 33e+09


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

    fast_num = 0
    for i,line in enumerate(lines):
        if "ITEM: ATOMS" in line:
            fast_num = i
            break
        
    for line in lines[fast_num+1:]:
        list = line.strip().split()
        x,y,z = map(float, list[2:5])
        atoms = np.array([x,y,z])
        return atoms

######## ファイル数の読み込み ########
with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/infilename.txt", 'r') as f_num:
    lines = f_num.readlines()
    file = lines[0].split()
    file_num = len(file)
    print(f"{file_num}")

######## supercell sizeの読み込み ########
disp_threshold = []
supercell_size = []
with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/supercell_size.txt", 'r') as f_size:
    lines = f_size.readlines()
    data_amount = len(lines)
    for i in range(data_amount-file_num, data_amount):
        size = np.array(lines[i].split(),dtype =float)
        supercell_size.append(size)
        half_size = size/2
        print(f'{half_size}')
        disp_threshold.append(half_size)

    
######## 座標データの読み込み ########
coordinate_of_defect = []
coordinate_of_perfect = []
coordinate_of_vacancy = []
atom_num = []
common_ids = []
for i in range(file_num):
    coordinate_p = read_perfect_data(f"/Users/kyou/Library/CloudStorage/Box-Box/dump/vacancy_Zr/perfect_lattice/perfect_lattice{i}.dump")
    coordinate_of_perfect.append(coordinate_p)
    coordinate_d = read_defect_data(f"/Users/kyou/Library/CloudStorage/Box-Box/dump/vacancy_Zr/defect_lattice/defect_lattice{i}.dump")
    coordinate_of_defect.append(coordinate_d)
    coordinate_v = read_vacancydata(f"/Users/kyou/Library/CloudStorage/Box-Box/dump/vacancy_Zr/deleted_atom/deleted_atom{i}.dump")
    coordinate_of_vacancy.append(coordinate_v)
    atom_num.append(len(coordinate_d))

    common_id = sorted(set(coordinate_p))
    common_ids.append(common_id)  

######## displacementの出力 ########
for i in range(file_num):
    with open(f"/Users/kyou/Library/CloudStorage/Box-Box/dump/vacancy_Zr/displacement/displacement{i}.inp", "w") as f_disp,\
        open(f"/Users/kyou/Library/CloudStorage/Box-Box/dump/vacancy_Zr/displacement/displacement_ovit{i}.xyz", "w") as f_disp_ovit:
        
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
                    print(f'{u[k]} {disp_threshold[i][k]}')

            f_disp_ovit.write(f"Zr {x[0]:.5f} {x[1]:.5f} {x[2]:.5f} {u[0]:.5f} {u[1]:.5f} {u[2]:.5f}\n")
            f_disp.write(f"{x[0]*conv_ang_to_m:.12e} {x[1]*conv_ang_to_m:.12e} {x[2]*conv_ang_to_m:.12e} {u[0]*conv_ang_to_m:.12e} {u[1]*conv_ang_to_m:.12e} {u[2]*conv_ang_to_m:.12e}\n")

            

