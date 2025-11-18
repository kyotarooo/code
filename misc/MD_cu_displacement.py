import numpy as np # type: ignore
import subprocess

#銅の材料定数の定義
E = 122.6e+09 # Pa 裳華房　材料力学　P.15
nu = 0.3 # ポアソン比　　裳華房　材料力学　P.15
g = E / (2 * (1 + nu)) # 剛性率
cutoff_radios = np.array(list(range(3, 31))) # 3から30までの整数（1刻み）
lattice_const = 3.615
force_dipole_factor = 1.0e-25

#ファイル名
lammpsscriptfilename = "lammps_cu.in"#ここで、lammpsのスクリプトファイルを指定している

with open(f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/{lammpsscriptfilename}",'w') as f:      
    f.write('log /Users/kyou/Library/CloudStorage/Box-Box/lammps/log.lammps\n')
    f.write('clear\n')#条件等全てリセット
    f.write('units metal\n')#単位系の指定　metalの場合は、Å
    f.write('dimension 3\n')
    f.write('boundary p p p\n')#周期境界条件
    f.write('atom_style atomic\n')
    f.write('lattice fcc 3.615\n')
    f.write('region my_region block 0 20 0 20 0 20\n')#my_regionは領域の名称、blochは直方体、単位は格子間隔
    f.write('create_box 1 my_region\n')#シミュレーションボックスの作成
    f.write('create_atoms 1 region my_region\n')#原子を生成
    f.write('mass 1 63.546\n')
    f.write('pair_style eam\n')#埋め込み原子法 Embedded Atom Method
    f.write('pair_coeff * * /Users/kyou/Library/CloudStorage/Box-Box/interatomic_potentials/Cu_u3.eam\n') 
    f.write('region vacancyspot sphere 10 10.1 10 0.5\n')#box中心に、半径０.１Åの球、格子間隔で定義
    f.write('group vacancyatom region vacancyspot\n')#group ID region R 事前の定義されたコマンドの領域内の原子をグループ化
    for n in cutoff_radios:
        f.write(f'region residual_applied_area_{n} sphere 10 10 10 {n/lattice_const}\n')
        f.write(f'group residual_applied_atom_{n} region residual_applied_area_{n}\n')
    f.write('dump 1 vacancyatom custom 1 /Users/kyou/Library/CloudStorage/Box-Box/lammps/deleted_atoms.dump id type x y z\n')
    f.write('dump_modify 1 format line "%d %d %.10f %.10f %.10f"\n')
    f.write('run 0\n')          # 0ステップ実行して、最初の状態を dump ファイルに書き出す
    f.write('undump 1\n')       # dump設定を解除
    f.write('delete_atoms group vacancyatom\n')
    f.write('dump 1 all custom 1 /Users/kyou/Library/CloudStorage/Box-Box/lammps/MD_Cu_perfect.dump id type x y z\n')
    f.write('run 0\n')          # 0ステップ実行して、最初の状態を dump ファイルに書き出す
    f.write('undump 1\n')       # dump設定を解除
    f.write('minimize 1e-8 1e-10 1000 10000\n')#minimize etol ftol maxiter maxeval  etol:stopping tolerance for energy (unitless),連続するステップのエネルギー差をエネルギーで割った値が、指定した許容値であるか否か
    #etol:内部応力の閾値
    f.write('dump 1 all custom 1 /Users/kyou/Library/CloudStorage/Box-Box/lammps/MD_Cu_defect.dump id type x y z\n')
    f.write('run 1000\n')#1000ステップで終了
    f.write('compute peratom all stress/atom NULL\n')
    for i in range(3, 31):
        f.write(f'compute sumstress{i} residual_applied_atom_{i} reduce sum '
                'c_peratom[1] c_peratom[2] c_peratom[3] c_peratom[4] c_peratom[5] c_peratom[6]\n')
    for i in range(3, 31):
        f.write(f'variable pxx{i} equal c_sumstress{i}[1]*{force_dipole_factor}\n')
        f.write(f'variable pyy{i} equal c_sumstress{i}[2]*{force_dipole_factor}\n')
        f.write(f'variable pzz{i} equal c_sumstress{i}[3]*{force_dipole_factor}\n')
        f.write(f'variable pxy{i} equal c_sumstress{i}[4]*{force_dipole_factor}\n')
        f.write(f'variable pxz{i} equal c_sumstress{i}[5]*{force_dipole_factor}\n')
        f.write(f'variable pyz{i} equal c_sumstress{i}[6]*{force_dipole_factor}\n')
        f.write(f'fix printstress{i} all print 1 "${{pxx{i}}} ${{pyy{i}}} ${{pzz{i}}} ${{pxy{i}}} ${{pxz{i}}} ${{pyz{i}}}" '
                f'file /Users/kyou/Library/CloudStorage/Box-Box/lammps/stress{i}.txt screen no\n')
    f.write('run 0\n')
    

       
subprocess.run(["/Users/kyou/venv/bin/lmp", "-in", f"/Users/kyou/Library/CloudStorage/Box-Box/lammps/{lammpsscriptfilename}"])#実行するコマンドをリストで指定し、ファイルのパスを指定して実行する
#subprocess.run: 外部のコマンドや他のプログラムを実行するための関数
#lmpは、LAMMPSの実行ファイルのパス、-inは、LAMMPSのスクリプトファイルを指定するオプション

def read_last_dump(filename):#dump: データとかを、そのままの形で出力すること
    with open(filename,'r') as f:
        lines=f.readlines()#readlines: ファイル全体を行ごとに分割したリストとして取得

    #最後のブロックのITEM: ATOMSを探す
    last_timestep_index=0

    for i,line in enumerate(lines):#enumarate:イテラブルオブジェクトの要素とインデックス番号を取得できる,このようにループ処理する際に用いる関数

        if "ITEM: ATOMS" in line:#if文の部分一致
            last_timestep_index=i
      
    atoms={}#dictの宣言。配列の宣言の場合、[]となる

    for atom_line in lines[last_timestep_index+1:]:#スライスは、[start:stop]では、start<= x < stopの範囲
        tokens = atom_line.strip().split()
        #strip():空白文字を削除する
        #split():()で分割する
        atom_id = int(tokens[0])#readlinesで読んだ情報は文字列だから、整数の型に変換
        x,y,z=map(float, tokens[2:5])#map:イテラブル関数の全要素を全て、指定した関数に適用できる
        atoms[atom_id]=np.array([x,y,z])
        #イテラブル：繰り返し可能なオブジェクト

    #atoms.sort(key=lambda x: x[0])#無名関数ラムダを用いている。lambda 引数：戻り値,慣習的によくXを用いる。
    #return np.array([pos for _, pos in atoms])#最初の要素は捨てて、二番目の要素を取り出す。
    return atoms


def read_fast_dump(filename):#dump: データとかを、そのままの形で出力すること
    with open(filename,'r') as f:
        lines=f.readlines()#readlines: ファイル全体を行ごとに分割したリストとして取得

    fast_timestep_index=100000

    for i,line in enumerate(lines):#enumarate:イテラブルオブジェクトの要素とインデックス番号を取得できる

        if "ITEM: ATOMS" in line:
            if i<fast_timestep_index:
                fast_timestep_index=i
      
    atoms={}#dictの宣言。配列の宣言の場合、[]となる

    for atom_line in lines[fast_timestep_index+1:]:#スライスは、[start:stop]では、start<= x < stopの範囲
        tokens = atom_line.strip().split()
        #strip():空白文字を削除する
        #split():()で分割する
        atom_id = int(tokens[0])#readlinesで読んだ情報は文字列だから、整数の型に変換
        x,y,z=map(float, tokens[2:5])#map:イテラブル関数の全要素を全て、指定した関数に適用できる
        atoms[atom_id]=np.array([x,y,z])
        #イテラブル：繰り返し可能なオブジェクト

        if "ITEM: ATOMS" in atom_line:
            break
    return atoms


def read_defect_coordinate(filename):
    with open(filename,'r') as f:
        lines = f.readlines()

    atoms = {}

    for i,line in enumerate(lines):

        if "ITEM: ATOMS" in line:
            for atom_line in lines[i+1:]:
                tokens = atom_line.strip().split()
                atom_id = int(tokens[0])#文字列から整数に変換
                x,y,z = map(float,tokens[2:5])
                atoms[atom_id] = np.array([x,y,z])
                break

    return atoms        

#データの読み込み
coordenate_of_defect = read_defect_coordinate("/Users/kyou/Library/CloudStorage/Box-Box/lammps/deleted_atoms.dump")
perfect = read_fast_dump("/Users/kyou/Library/CloudStorage/Box-Box/lammps/MD_Cu_perfect.dump")
number_of_atom_perfect=len(perfect)
defect = read_last_dump("/Users/kyou/Library/CloudStorage/Box-Box/lammps/MD_Cu_defect.dump")
number_of_atom_defect=len(defect)

#set:集合を作成する関数、setは重複を許さない(keyの重複を許さない)、集合の演算ができる
#sorted:ソートする関数、setをソートしてリストに変換
common_ids = sorted(set(perfect))
deleted_ids = list(coordenate_of_defect.keys())

print(f'{number_of_atom_perfect},{number_of_atom_defect}')
print(f'{coordenate_of_defect}')



scale_factor=100
with open("/Users/kyou/Library/CloudStorage/Box-Box/lammps/MD_Cu_displacement.xyz","w") as f1:#ovit用
    with open("/Users/kyou/Library/CloudStorage/Box-Box/code/coordinate/data.inp","w") as f2:#green関数用
        f1.write(f"{number_of_atom_defect}\n")
        # 拡張XYZのヘッダ行（各列の意味を定義）
        # species: 元素記号, pos: 座標, disp: 変位ベクトル（これがOVITOでベクトル表示される）
        f1.write("Properties=species:S:1:pos:R:3:disp:R:3\n")

        unit_factor = 1.0e-10#単位はメートル

        f2.write(f"{g} {nu}\n")
        xc1,xc2,xc3 = coordenate_of_defect[deleted_ids[0]]
        f2.write(f"{xc1*unit_factor} {xc2*unit_factor} {xc3*unit_factor}\n")
        f2.write(f"{number_of_atom_defect}\n")
    
        for atom_id in common_ids:
            x,y,z=perfect[atom_id]
            u=defect[atom_id]-perfect[atom_id]
            ux,uy,uz=u
            scale_factor=100
            f1.write(f"Cu {x:.5f} {y:.5f} {z:.5f} {ux*scale_factor:.5f} {uy*scale_factor:.5f} {uz*scale_factor:.5f}\n")#f文字列の桁数指
            scale_factor=1
            f2.write(f"{x*unit_factor:.12e} {y*unit_factor:.12e} {z*unit_factor:.12e} {ux*unit_factor:.12e} {uy*unit_factor:.12e} {uz*unit_factor:.12e}\n")
            #少数点以下の桁数指定



