#!/bin/zsh

######## 出力ファイルの入力 ########
echo -n 'material: '
read material  
echo -n 'potential: '
read potential
echo -n 'defect: '
read defect
echo -n 'method: '
read method

####### 出力ファイルのパス ########
output_base_dir=~/Library/CloudStorage/Box-Box/output
output_dir=$material/$defect/$potential/$method
code_dir=~/Library/CloudStorage/Box-Box/code/forcedipole_${material}_${defect}/$potential/$method
full_output_dir=$output_base_dir/$output_dir
echo "-------------------------------------------------"
echo "Output Path: $full_output_dir"
echo "-------------------------------------------------"
echo "Press [Enter]"
read  

####### make directory ########
# base output file
rm -rf $full_output_dir/fullpath.txt
mkdir -p "$full_output_dir"
echo "$full_output_dir" > ${full_output_dir}/fullpath.txt
# dump file (ovit)
rm -rf $full_output_dir/dump
mkdir -p "$full_output_dir/dump"
# lammps log file
rm -rf $full_output_dir/log_file
mkdir -p "$full_output_dir/log_file"

# supercell 
rm -f $full_output_dir/supercell.txt

# perfect lattice coordinate
for i in {0..11}; do
    rm -rf "$full_output_dir/dump/perfect_lattice/perfect_lattice_$i"
    mkdir -p "$full_output_dir/dump/perfect_lattice/perfect_lattice_$i"
done
  
# defect lattice coordinate
for i in {0..11}; do
    rm -rf "$full_output_dir/dump/defect_lattice/defect_lattice_$i"
    mkdir -p "$full_output_dir/dump/defect_lattice/defect_lattice_$i"
done

# deleted atom coordinate
for i in {0..11}; do
    rm -rf "$full_output_dir/dump/deleted_atom/deleted_atom_$i"
    mkdir -p "$full_output_dir/dump/deleted_atom/deleted_atom_$i"
done

# displacement data
for i in {0..11}; do
    rm -rf "$full_output_dir/dump/displacement/displacement_$i"
    mkdir -p "$full_output_dir/dump/displacement/displacement_$i"
done

# displacement ovit_data
for i in {0..11}; do
    rm -rf "$full_output_dir/dump/displacement_ovit/displacement_ovit_$i"
    mkdir -p "$full_output_dir/dump/displacement_ovit/displacement_ovit_$i"
done

# force dipole file
for i in {0..11}; do
    rm -rf "$full_output_dir/dump/force_dipole/force_dipole_$i"
    mkdir -p "$full_output_dir/dump/force_dipole/force_dipole_$i"
done

# 4hsic_vacancy(0 - 11)s
for i in {0..11}; do
    rm -rf "$full_output_dir/4hsic_vacancy/atom_type_delete_$i"
    mkdir -p "$full_output_dir/4hsic_vacancy/atom_type_delete_$i"
done

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

export ATOM=0 #どの位置の原子を削除するか ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←


######## make sic structure ########

echo "-------------------------------------------------"
echo "make 4h-sic structure"
echo "-------------------------------------------------"
echo "Press [Enter]"
read

python3 "$code_dir/sic.py"

######## delete atom ########

echo "-------------------------------------------------"
echo "delete 1 atom"
echo "-------------------------------------------------"
echo "Press [Enter]"
read

python3 "$code_dir/new_delete_atom.py"

######## run rammps script ########
# lmpファイル
LAMMPS=/Users/kyou/Library/CloudStorage/Box-Box/lammps-29Oct20/src/lmp_mpi

# lammps script file
lammps_script=minimize.lammps
echo "-------------------------------------------------"
echo "run lammps script"
echo "-------------------------------------------------"
echo "Press [Enter]"
read

mpirun -np 8 "$LAMMPS" -var output_dir $OUTPUT_PATH -in "$code_dir/$lammps_script" 


######## read data ########
echo "-------------------------------------------------"
echo "read data"
echo "-------------------------------------------------"
echo "Press [Enter]"
read
python3 "$code_dir/readdata.py"

######## green function ########
echo "-------------------------------------------------"
echo "calculate green function"
echo "-------------------------------------------------"
echo "Press [Enter]"
read
python3 "$code_dir/green.py"

######## visualization ########
echo "-------------------------------------------------"
echo "plot"
echo "-------------------------------------------------"
echo "Press [Enter]"
read
python3 "$code_dir/visualization_forcedipole.py"

