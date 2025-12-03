#!/bin/zsh 
#PBS -l nodes=intal01:ppn=7
#PBS -e stderr.txt 
#PBS -o stdout.txt 
#PBS -N p3

cd $PBS_O_WORKDIR 
NPROCS=`wc -l < $PBS_NODEFILE` 
cd `pwd` 

export LD_LIBRARY_PATH="/home/common/intel/oneapi/tbb/2021.9.0/env/../lib/intel64/gcc4.8:/home/common/intel/oneapi/mpi/2021.9.0//libfabric/lib:/home/common/intel/oneapi/mpi/2021.9.0//lib/release:/home/common/intel/oneapi/mpi/2021.9.0//lib:/home/common/intel/oneapi/mkl/2023.1.0/lib/intel64:/home/common/intel/oneapi/itac/2021.9.0/slib:/home/common/intel/oneapi/ipp/2021.8.0/lib/intel64:/home/common/intel/oneapi/ippcp/2021.7.0/lib/intel64:/home/common/intel/oneapi/ipp/2021.8.0/lib/intel64:/home/common/intel/oneapi/dnnl/2023.1.0/cpu_dpcpp_gpu_dpcpp/lib:/home/common/intel/oneapi/debugger/2023.1.0/gdb/intel64/lib:/home/common/intel/oneapi/debugger/2023.1.0/libipt/intel64/lib:/home/common/intel/oneapi/debugger/2023.1.0/dep/lib:/home/common/intel/oneapi/dal/2023.1.0/lib/intel64:/home/common/intel/oneapi/compiler/2023.1.0/linux/lib:/home/common/intel/oneapi/compiler/2023.1.0/linux/lib/x64:/home/common/intel/oneapi/compiler/2023.1.0/linux/lib/oclfpga/host/linux64/lib:/home/common/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin:/home/common/intel/oneapi/ccl/2021.9.0/lib/cpu_gpu_dpcpp" 

DD_PATH="/home/kino/ElementBased-orowantest/ebpdd" 
log_name="dd.log" 

${DD_PATH} ./ > ${log_name} 
