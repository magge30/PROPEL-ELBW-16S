#!/bin/bash -l
#SBATCH -A <SNICProject>
#SBATCH -p core
#SBATCH -n 20
#SBATCH -J bbduk
#SBATCH -t 5:00:00
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

# make a loop to trim samples with bbduk
ls *fastq.gz | cut -f1-2 -d "_" > samples

for sample in $(cat samples); \
do echo "On sample: $sample"; bbduk.sh in="$sample"_L001_R1_001.fastq.gz in2="$sample"_L001_R2_001.fastq.gz \
out="$sample"_R1_trimmed1.fastq.gz out2="$sample"_R2_trimmed1.fastq.gz \
literal=CCTACGGGNGGCWGCAG,GACTACHVGGGTATCTAATCC k=17 ordered=t \
ktrim=l rcomp=f mm=f copyundefined ftm=5 qtrim=l trimq=35 minlen=100; \
done 2> trimming1_stats.txt

# make a loop to trim samples with bbduk
ls *trimmed1.fastq.gz | cut -f1-2 -d "_" > trimmed1_samples

for sample in $(cat trimmed1_samples); \
do echo "On sample: $sample"; bbduk.sh in="$sample"_R1_trimmed1.fastq.gz in2="$sample"_R2_trimmed1.fastq.gz \
out="$sample"_R1_trimmed.fastq.gz out2="$sample"_R2_trimmed.fastq.gz \
literal=CCTACGGGNGGCWGCAG,GACTACHVGGGTATCTAATCC k=17 ordered=t \
ktrim=l rcomp=f mm=f copyundefined ftm=5 qtrim=rl trimq=30 minlen=100; \
done 2> trimming_stats.txt

#quality check trimmed
mkdir quality_check_trimmed
fastqc -o quality_check_trimmed *trimmed.fastq.gz
cd quality_check_trimmed
multiqc .
