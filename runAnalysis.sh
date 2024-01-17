#!/bin/bash -l
# How long should I job run for
#SBATCH --time=48:00:00
# Number of CPU cores, in this case 1 core
#SBATCH --ntasks=40
# Number of compute nodes to use
#SBATCH --nodes=1
# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Run the analysis in the singularity container - binding the current working directory
export OMP_NUM_THREADS=40
module load singularity
singularity run -B `pwd -P` runAnalysis.img
