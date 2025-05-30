#!/bin/bash

# Base script generated by NERSC Batch Script Generator on https://iris.nersc.gov/jobscript.html

#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J ap_train
#SBATCH --mail-user=lschlueter@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -A m2676
#SBATCH -t 0:10:00

# Default values if no arguments are provided
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application: 
export LEGEND_DATA_CONFIG=/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json
srun -n 1 -c 256 --cpu_bind=cores /global/common/software/nersc/n9/julia/1.10.4/bin/julia --project=/global/homes/l/lschl/.julia/environments/asic-dev/ ./02_AP_train.jl
