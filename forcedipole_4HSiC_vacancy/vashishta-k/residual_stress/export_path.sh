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
output_dir_R=$material/$defect/$potential/residual_stress
code_dir=~/Library/CloudStorage/Box-Box/code/forcedipole_${material}_${defect}/$potential/$method
full_output_dir=$output_base_dir/$output_dir

######## 削除する原子の指定 ######## ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
target_atom=0

####### 材料定数の定義 ########
export LATTICE_CONST=3.076599647 # ang
export LATTICE_TYPE="Diamond"

export OUTPUT_PATH="$full_output_dir"
export DIPOLE_R_PATH="$output_base_dir/$output_dir_R"

export ATOM=$target_atom 

echo "-------------------------------------------------"
echo "target atom: $target_atom"
echo "-------------------------------------------------"
echo "Press [Enter]"
read  


