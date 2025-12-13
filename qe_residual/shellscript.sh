#!/bin/zsh
# =========================================================================================
# 注意
# burgersの変更
# 削除する原子が共通になるように
# =========================================================================================

#!/bin/zsh

######## 出力ファイルの入力 ########
echo -n 'material: '
read material 
echo -n 'defect: '
read defect
echo -n 'method: '
read method

####### 出力ファイルのパス ########
output_base_dir=~/Library/CloudStorage/Box-Box/output
output_dir=$method/$material/$defect
code_dir=~/Library/CloudStorage/Box-Box/code/$method
full_output_dir=$output_base_dir/$output_dir
echo "-------------------------------------------------"
echo "Output Path: $full_output_dir"
echo "-------------------------------------------------"
echo "Press [Enter]"
read  

######## 削除する原子の指定 ######## ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←minimize.pyも変更する
target_atom=5

####### make directory ########
# base output file
rm -rf $full_output_dir/fullpath.txt
mkdir -p "$full_output_dir"
echo "$full_output_dir" > ${full_output_dir}/fullpath.txt

# lammps log file
rm -rf $full_output_dir/log_file
mkdir -p "$full_output_dir/log_file"

# supercell 
rm -f $full_output_dir/supercell.txt


########## ディレクトリ作成・削除 ##########

# force dipole file
rm -rf "$full_output_dir/force_dipole/force_dipole_${target_atom}"
mkdir -p "$full_output_dir/force_dipole/force_dipole_${target_atom}"

# perfect lattice coordinate
rm -rf "$full_output_dir/dump/perfect_lattice/perfect_lattice_${target_atom}"
mkdir -p "$full_output_dir/dump/perfect_lattice/perfect_lattice_${target_atom}"

# defect lattice coordinate
rm -rf "$full_output_dir/dump/defect_lattice/defect_lattice_${target_atom}"
mkdir -p "$full_output_dir/dump/defect_lattice/defect_lattice_${target_atom}"

# deleted atom coordinate
rm -rf "$full_output_dir/dump/deleted_atom/deleted_atom_${target_atom}"
mkdir -p "$full_output_dir/dump/deleted_atom/deleted_atom_${target_atom}"

# # displacement data
# rm -rf "$full_output_dir/dump/displacement/displacement_${target_atom}"
# mkdir -p "$full_output_dir/dump/displacement/displacement_${target_atom}"

# # displacement (ovit)
# rm -rf "$full_output_dir/dump/displacement_ovit/displacement_ovit_${target_atom}"
# mkdir -p "$full_output_dir/dump/displacement_ovit/displacement_ovit_${target_atom}"

# 4hsic_vacancy
rm -rf "$full_output_dir/4hsic_vacancy/atom_type_delete_${target_atom}"
mkdir -p "$full_output_dir/4hsic_vacancy/atom_type_delete_${target_atom}"


# 4hsic
rm -rf $full_output_dir/4hsic
mkdir -p $full_output_dir/4hsic
# 4hsic_q
rm -rf $full_output_dir/4hsic_q
mkdir -p $full_output_dir/4hsic_q
# include
rm -rf $full_output_dir/include
mkdir -p $full_output_dir/include


# delite same files
rm -rf $full_output_dir/N_atom_perfect.txt
rm -rf $full_output_dir/N_atom_perfect.txt

####### 材料定数の定義 ########
export LATTICE_CONST=3.076599647 # ang
export LATTICE_TYPE="Diamond"

export OUTPUT_PATH="$full_output_dir"
export DIPOLE_R_PATH="$code_dir/residual_stress/force_dipole"

export ATOM=$target_atom #どの位置の原子を削除するか ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
export Burgers=3.0e-10 #m ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←


######## make sic structure ########

echo "-------------------------------------------------"
echo "make 4h-sic structure"
echo "-------------------------------------------------"
echo "Press [Enter]"
read

python3 "$code_dir/sic.py"

######## make Vsi-4H ########

echo "-------------------------------------------------"
echo "make Vsi-4H"
echo "-------------------------------------------------"
echo "Press [Enter]"
read

python3 "$code_dir/vsi4hgen.py"


