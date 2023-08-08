#!/bin/bash

#SBATCH --nodes=15
#SBATCH --time=12:00:00
#SBATCH --account=hfm
#SBATCH --job-name=crz_RANS_benchmark_ADM
#SBATCH --output=out.%x_%j
#SBATCH --mail-user=shen.huang@nrel.gov
#SBATCH --mail-type=BEGIN,END,FAIL
export SPACK_MANAGER=${HOME}/exawind_spack/spack-manager
source ${SPACK_MANAGER}/start.sh
spack-start
spack env activate -d ${HOME}/exawind_spack/spack-manager/environments/myenv
spack load amr-wind


srun -n 540 amr_wind crz_RANS_benchmark.inp &> log
