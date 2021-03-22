#!/bin/bash -l
#SBATCH -A <SNICproject>
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J quality_check_raw
#SBATCH -t 03:00:00
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>

# load modules
module load bioinfo-tools
module load FastQC/0.11.5
module load MultiQC/1.7

# go to working directory
cd /your/working/directory

# run fastqc and multiqc
mkdir quality_check_raw
fastqc -o quality_check_raw *.fastq.gz
cd quality_check_raw
multiqc .
