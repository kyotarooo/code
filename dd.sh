#!/bin/zsh
#PBS -N AmdRa01
#PBS -l nodes=amdve06:ppn=1
#PBS -e stderr.txt
#PBS -o stdout.txt

cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
# cd `pwd`

log_name="dd.log"

DD_PATH="/home/kyou/InclusionStress-20251124/incstrgen"
${DD_PATH} ./ > ${log_name}
