#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=sharpen
#SBATCH --time=0:01:0
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=1
# You should normally set this to your own budget, e.g. m25oc-s1234567
#SBATCH --account=m25oc
#SBATCH --partition=standard
#SBATCH --qos=short

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

# Launch the parallel job
#   Using 4 MPI processes on 1 node
#   srun picks up the distribution from the sbatch options
srun --unbuffered --distribution=block:block --hint=nomultithread ./sharpen

