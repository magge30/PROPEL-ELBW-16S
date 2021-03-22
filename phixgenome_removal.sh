#!/bin/bash -l
#SBATCH -A <SNICProject>
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J phix_removal
#SBATCH -t 3:00:00
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>

# load modules
module load bioinfo-tools
module load FastQC/0.11.5
module load MultiQC/1.7
module load bbmap/38.08

# go to working directory
cd /your/working/directory

# make a loop to reads mapping to the phix genome with bbduck
ls *trimmed.fastq.gz | cut -f1-2 -d "_" > trimmed_samples

for sample in $(cat trimmed_samples); \
do echo "On sample: $sample"; bbduk.sh in="$sample"_R1_trimmed.fastq.gz in2="$sample"_R2_trimmed.fastq.gz \
out="$sample"_R1_trimmed_phixunmatch.fastq.gz out2="$sample"_R2_trimmed_phixunmatch.fastq.gz \
outm="$sample"_R1_trimmed_phixmatch.fastq.gz outm2="$sample"_R2_trimmed_phixmatch.fastq.gz \
ref=genome.fa k=31 hdist=1; \
done 2> phix_stats.txt

#Quality check trimmed_unmatch
mkdir quality_check_phixunmatch
fastqc -o quality_check_phixunmatch *trimmed_phixunmatch.fastq.gz
cd quality_check_phixunmatch
multiqc .
