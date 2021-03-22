#!/bin/bash -l
#SBATCH -A <SNICProject>
#SBATCH -p node
#SBATCH -N 4
#SBATCH -J dada2_silva
#SBATCH -t 72:00:00
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>

# load modules
module load R
module load R_packages/3.5.2
module load gcc/8.2.0
module load openmpi/3.1.1

# go to working directory where Rscript is located
cd /your/working/directory

export PATH=$PATH:/home/magali/R/x86_64-pc-linux-gnu-library/3.5/snow
mpirun -np 80 RMPISNOW --save < dada2_rscript.R

