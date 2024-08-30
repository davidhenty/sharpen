#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=sharpen
#SBATCH --time=0:01:0
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --account=tc063
#SBATCH --partition=standard
#SBATCH --qos=short

# Ensure cpus-per-task is exported to srun
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Load the default HPE MPI environment
module load mpt
module load intel-20.4/compilers

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

# Set the number of threads to the CPUs per task 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Launch the parallel job
#   Using 36 threads per node
#   srun picks up the distribution from the sbatch options
srun --cpu-bind=cores ./sharpen

