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
code_dir=~/Library/CloudStorage/Box-Box/code/forcedepole_${material}_${defect}/$potential/$method
full_output_dir=$output_base_dir/$output_dir
echo "-------------------------------------------------"
echo "Output Path: $full_output_dir"
echo "-------------------------------------------------"
echo "Press [Enter]"
read  

####### make directory ########
# base output file
mkdir -p "$full_output_dir"
echo "$full_output_dir" > ${full_output_dir}/fullpath.txt
# lammps .in file
mkdir -p "$full_output_dir/in_file"
# lammps log file
mkdir -p "$full_output_dir/log_file"
# force dipole file
mkdir -p $full_output_dir/force_dipole
# sumstress file
rm -f $full_output_dir/N_atom_defect.txt
rm -f $full_output_dir/N_atom_perfect.txt
rm -f $full_output_dir/sumstress
mkdir -p $full_output_dir/sumstress
# dump file (ovit)
#mkdir -p "$full_output_dir/dump"
# deleted atom coordinate
#mkdir -p $full_output_dir/dump/deleted_atom
# perfect lattice coordinete
#mkdir -p $full_output_dir/dump/perfect_lattices
# defect lattice coordinate
#mkdir -p $full_output_dir/dump/defect_lattice
# displacement data
#mkdir -p $full_output_dir/dump/displacement
# displacement ovit data
#mkdir -p $full_output_dir/displacement_ovit

####### make lammps script ########
make_infile_script=make_${material}_${defect}.py
echo "-------------------------------------------------"
echo "Make lammps script: $code_dir/$make_infile_script"
echo "-------------------------------------------------"
echo "Press [Enter]"
read
export OUTPUT_PATH="$full_output_dir"
python3 "$code_dir/$make_infile_script"

######## run rammps script ########
# lmpファイル
LAMMPS=~/lammps_install/bin/lmp
# lammps scriptフォルダ
INPUT_DIR=$full_output_dir/in_file
# lammps script name
LIST=$full_output_dir/infilename.txt
echo "-------------------------------------------------"
echo "run lammps script"
echo "-------------------------------------------------"
echo "Press [Enter]"
read
# リストを1行ずつ読み取って実行
for lammps_script in $(cat "$LIST"); do
    echo "===== Running $lammps_script ====="
    mpirun -np 1 "$LAMMPS" -in "$INPUT_DIR/$lammps_script"
done

######## subtraction ########
echo "-------------------------------------------------"
echo "subtraction"
echo "-------------------------------------------------"
echo "Press [Enter]"
read
python3 "$code_dir/subtraction.py"


######## visualization ########
echo "-------------------------------------------------"
echo "plot"
echo "-------------------------------------------------"
echo "Press [Enter]"
read
python3 "$code_dir/visualization_forcedipole.py"




