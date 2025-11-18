import numpy as np # type: ignore
import os

######## 環境変数の取得 ########
output_dir = os.environ.get("OUTPUT_PATH")
if output_dir is None:
    raise ValueError("環境変数 OUTPUT_PATH error")

######### 材料定数の取得 ############
lattice_const = os.environ.get("LATTICE_CONST") # ang
lattice_const = float(lattice_const)
lattice_type = os.environ.get("LATTICE_TYPE")

######### 単位変換係数 ########
conv_eV_to_J = 1.602176634E-19
conv_bars_to_Pa = 1.0E+5
conv_Vang_to_Vm = 1.0E-30
conv_barVang_to_NVm = conv_bars_to_Pa * conv_Vang_to_Vm
conv_J_to_eV = 1/conv_eV_to_J

######## supercell のサイズ ########
units_cell_x = [2, 4, 6, 8, 10, 12]
units_cell_y = [2, 4, 6, 8, 10, 12]
units_cell_z = [2, 4, 6, 8, 10, 12]
atom_num = []
filename = []

######## ファイルの名前 ########
for i in range(len(units_cell_x)):
    num = f"{units_cell_x[i]}x{units_cell_y[i]}x{units_cell_z[i]}"
    name = f"lammps_{num}.in"
    atom_num.append(num)
    filename.append(name)


######## ファイル名の出力 ########  
with open(f"{output_dir}/infilename.txt", 'w') as f_num:
    for i in range(len(units_cell_x)):
        f_num.write(f"{filename[i]} ")

######## lammps script #########
for i in range(len(units_cell_x)):
    with open(f"{output_dir}/in_file/{filename[i]}", 'w') as f:  
        f.write(f'log {output_dir}/log_file/log.lammps\n')
        f.write('clear\n')
        f.write('units metal\n')
        f.write('boundary p p p\n')
        f.write('atom_style atomic\n')

        # HCP格子
        f.write(f'lattice custom {lattice_const} a1 1.0 0.0 0.0 a2 -0.5 0.8660254 0.0 a3 0.0 0.0 1.598 '
                'basis 0.0 0.0 0.0 basis 0.6666667 0.3333333 0.5\n')#c/a ratio 
        
        # supercellの作成
        f.write(f'region supercell{i} block 0 {units_cell_x[i]} 0 {units_cell_y[i]} 0 {units_cell_z[i]}\n')
        f.write(f'create_box 1 supercell{i}\n')
        f.write(f'create_atoms 1 region supercell{i}\n')
        #f.write(f'mass 1 {mass}\n')

        #原子数の確認
        f.write(f'variable N equal count(all)\n')
        f.write(f'print "${{N}}" append {output_dir}/N_atom_perfect.txt screen no\n')

        #lattce sizeの出力
        f.write(f'variable ax equal xlat\n')
        f.write(f'variable ay equal ylat\n')
        f.write(f'variable az equal zlat\n')
        f.write(f'print "${{ax}} ${{ay}} ${{az}}" append {output_dir}/lattice_size.txt screen no\n')

        #supercell sizeの出力
        f.write(f'variable lx equal lx\n')
        f.write(f'variable ly equal ly\n')
        f.write(f'variable lz equal lz\n')
        f.write(f'print "${{lx}} ${{ly}} ${{lz}}" append {output_dir}/supercell_size.txt screen no\n')
        
        # ポテンシャル
        f.write('pair_style eam/alloy\n')
        f.write('pair_coeff * * /Users/kyou/Library/CloudStorage/Box-Box/interatomic_potentials/Zr_3.eam.alloy Zr\n')

        # #緩和
        # f.write("min_style cg\n")
        # f.write("minimize 1.0e-10 1.0e-10 10000 100000\n")
        # f.write("run 1\n")

        # 応力の計算
        # f.write('compute pe_peratom all stress/atom NULL\n')
        # f.write(f'compute pe_sumstress{i} all reduce sum c_pe_peratom[1] c_pe_peratom[2] c_pe_peratom[3] c_pe_peratom[4] c_pe_peratom[5] c_pe_peratom[6]\n')
        # f.write('run 1\n')
        #*****
        #stress/atomで、原子ごとのローカルな応力を求め、reduce sumで合計をとることによって、グローバルな値にしている
        #*****

        # 応力の出力(perfect lattice (J))
        # f.write(f"variable V equal vol\n")
        # f.write(f"variable pe_pxx_J equal c_pe_sumstress{i}[1]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable pe_pyy_J equal c_pe_sumstress{i}[2]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable pe_pzz_J equal c_pe_sumstress{i}[3]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable pe_pxy_J equal c_pe_sumstress{i}[4]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable pe_pxz_J equal c_pe_sumstress{i}[5]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable pe_pyz_J equal c_pe_sumstress{i}[6]*{conv_barVang_to_NVm}\n")
        # f.write(f'print "{i} ${{pe_pxx_J}} ${{pe_pyy_J}} ${{pe_pzz_J}} ${{pe_pxy_J}} ${{pe_pxz_J}} ${{pe_pyz_J}}" append /Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_J/perfect_J_sumstress.txt screen no\n')
        # f.write('run 0\n')

         # 応力の出力(perfect lattice (eV))
        # f.write(f"variable V equal vol\n")
        # f.write(f"variable pe_pxx_eV equal c_pe_sumstress{i}[1]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable pe_pyy_eV equal c_pe_sumstress{i}[2]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable pe_pzz_eV equal c_pe_sumstress{i}[3]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable pe_pxy_eV equal c_pe_sumstress{i}[4]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable pe_pxz_eV equal c_pe_sumstress{i}[5]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable pe_pyz_eV equal c_pe_sumstress{i}[6]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f'print "{i} ${{pe_pxx_eV}} ${{pe_pyy_eV}} ${{pe_pzz_eV}} ${{pe_pxy_eV}} ${{pe_pxz_eV}} ${{pe_pyz_eV}}" append /Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_eV/perfect_eV_sumstress.txt screen no\n')
        # f.write('run 0\n')

        # 空孔の作成
        vacancy_coordinate = [units_cell_x[i]/2, units_cell_x[i]/2 + 0.1, units_cell_z[i]/2]
        f.write(f'region vacancy_area{i} sphere {vacancy_coordinate[0]} {vacancy_coordinate[1]} {vacancy_coordinate[2]} 0.37\n')
        f.write(f'group vacancy_atom{i} region vacancy_area{i}\n')

        # dump (deleted atom)
        f.write(f'dump deleted_atoms vacancy_atom{i} custom 1 {output_dir}/dump/deleted_atom/deleted_atom{i}.dump id type x y z\n')
        f.write('dump_modify deleted_atoms format line "%d %d %.15e %.15e %.15e"\n')
        f.write('run 0\n')

        # 空孔の削除
        f.write(f'delete_atoms group vacancy_atom{i}\n')
        f.write(f'variable N equal count(all)\n')
        f.write(f'print "${{N}}" append {output_dir}/N_atom_defect.txt screen no\n')
        f.write(f'run 0\n')

        # dump (perfect lattice)
        f.write(f'dump perfect all custom 1 {output_dir}/dump/perfect_lattice/perfect_lattice{i}.dump id type x y z\n') 
        f.write('dump_modify perfect format line "%d %d %.15e %.15e %.15e"\n')
        f.write("run 0\n")

        # 緩和
        f.write("min_style cg\n")
        f.write("minimize 1.0e-10 1.0e-10 10000 100000\n")
        f.write("run 1\n") 

        # dump (defect lattice)
        f.write(f'dump defect all custom 1 {output_dir}/dump/defect_lattice/defect_lattice{i}.dump id type x y z\n') 
        f.write('dump_modify defect format line "%d %d %.15e %.15e %.15e"\n')
        f.write("run 0\n")

        #ポテンシャルの出力
        #f.write("compute peratom all pe/atom\n")
        #f.write(f'dump 1 all custom 1 /Users/kyou/Library/CloudStorage/Box-Box/lammps/defect{i}.dump id type x y z c_peratom\n')
        
        # 応力の計算
        # f.write('compute de_peratom all stress/atom NULL\n')
        # f.write(f'compute de_sumstress{i} all reduce sum c_de_peratom[1] c_de_peratom[2] c_de_peratom[3] c_de_peratom[4] c_de_peratom[5] c_de_peratom[6]\n')
        # f.write('run 1\n')
      
        # 応力の出力(defect lattice (J))
        # f.write(f"variable de_pxx_J equal c_de_sumstress{i}[1]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable de_pyy_J equal c_de_sumstress{i}[2]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable de_pzz_J equal c_de_sumstress{i}[3]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable de_pxy_J equal c_de_sumstress{i}[4]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable de_pxz_J equal c_de_sumstress{i}[5]*{conv_barVang_to_NVm}\n")
        # f.write(f"variable de_pyz_J equal c_de_sumstress{i}[6]*{conv_barVang_to_NVm}\n")
        # f.write(f'print "{i} ${{de_pxx_J}} ${{de_pyy_J}} ${{de_pzz_J}} ${{de_pxy_J}} ${{de_pxz_J}} ${{de_pyz_J}}" append /Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_J/defect_J_sumstress.txt screen no\n')
        # f.write('run 0\n')

        # 応力の出力(defect lattice (eV))
        # f.write(f"variable de_pxx_eV equal c_de_sumstress{i}[1]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable de_pyy_eV equal c_de_sumstress{i}[2]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable de_pzz_eV equal c_de_sumstress{i}[3]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable de_pxy_eV equal c_de_sumstress{i}[4]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable de_pxz_eV equal c_de_sumstress{i}[5]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f"variable de_pyz_eV equal c_de_sumstress{i}[6]*{conv_barVang_to_NVm}*{conv_J_to_eV}\n")
        # f.write(f'print "{i} ${{de_pxx_eV}} ${{de_pyy_eV}} ${{de_pzz_eV}} ${{de_pxy_eV}} ${{de_pxz_eV}} ${{de_pyz_eV}}" append /Users/kyou/Library/CloudStorage/Box-Box/lammps/pointdefects_Zr/force_dipole_eV/defect_eV_sumstress.txt screen no\n')
        # f.write('run 0\n')






        
        