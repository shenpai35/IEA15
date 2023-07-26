#!/bin/bash

#SBATCH --nodes=10
#SBATCH --time=15:00:00
#SBATCH --account=hfm
#SBATCH --job-name=ws08_0_yaw_joukowsky_disk
#SBATCH --output=out.%x_%j
#SBATCH --mail-user=shen.huang@nrel.gov
#SBATCH --mail-type=BEGIN,END,FAIL
export SPACK_MANAGER=${HOME}/exawind_spack/spack-manager
source ${SPACK_MANAGER}/start.sh
spack-start
spack env activate -d ${HOME}/exawind_spack/spack-manager/environments/myenv
spack load amr-wind


srun -n 360 amr_wind 0_yaw_joukowsky_disk.inp &> log
