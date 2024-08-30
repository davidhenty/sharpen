#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=sharpen
#SBATCH --time=0:10:0
#SBATCH --ntasks=1
#SBATCH --account=tc063
#SBATCH --partition=standard
#SBATCH --qos=short

# Load the Python evirnment
module load python/3.9.13

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
 
# Launch the serial job
#   Using 1 thread
srun --cpu-bind=cores python ./sharpen.py

