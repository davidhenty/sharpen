#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=sharpen
#SBATCH --time=0:01:0
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=standard
#SBATCH --qos=short

# Ensure cpus-per-task is exported to srun
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Set the number of threads to the CPUs per task 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Ensure threads are assigned to cores appropriately
export OMP_PLACES=cores

# Launch the parallel job
#   Using 32 threads per node
#   srun picks up the distribution from the sbatch options
srun --unbuffered --distribution=block:block --hint=nomultithread ./sharpen

