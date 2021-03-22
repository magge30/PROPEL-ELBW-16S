# PROPEL-ELBW-16
Scripts to produce the analysis and figures for PROPEL 16S paper

Worklfow from raw sequences to ASV table: <br> 
#1.raw sequences quality check (fastqc and multiqc for raw reads)
quality_check_raw.sh

#2.quality trimming and filtering + quality check afterwards
bbduk_trimming.sh #qtrim=l trimq=35 + qtrim=rl trimq=30 minlen=100 + quality check

#3.remove phix genome NC_001422.1 from trimmed sequences + quality cehck afterwards
phixgenome_removal.sh

#4.denoising and taxonomic assigment (pseudo pooling, silva database, minBoot=80)
Trimmed sequences are processed following the DADA2 Workflow with dada2 version 1.10.1 (https://benjjneb.github.io/dada2/tutorial.html)
dada2_rscript.R
dada2.sh

Worklfow for statistical analysis:


