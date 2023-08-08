#!/bin/bash

#SBATCH --nodes=10
#SBATCH --time=12:00:00
#SBATCH --account=hfm
#SBATCH --job-name=crz_benchmark_lawrence_ADM
#SBATCH --output=out.%x_%j
#SBATCH --mail-user=shen.huang@nrel.gov
#SBATCH --mail-type=BEGIN,END,FAIL
export SPACK_MANAGER=${HOME}/exawind_spack/spack-manager
source ${SPACK_MANAGER}/start.sh
spack-start
spack env activate -d ${HOME}/exawind_spack/spack-manager/environments/myenv
spack load amr-wind


srun -n 360 amr_wind crz_benchmark_lawrence.inp &> log
