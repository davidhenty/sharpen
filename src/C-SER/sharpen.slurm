#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=sharpen
#SBATCH --time=0:01:0
#SBATCH --ntasks=1
#SBATCH --account=tc063
#SBATCH --partition=standard
#SBATCH --qos=short

# Load the default HPE MPI environment
module load mpt
module load intel-20.4/compilers

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
 
# Launch the serial job
#   Using 1 thread
srun --cpu-bind=cores ./sharpen

