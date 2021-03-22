# PROPEL-ELBW-16
Scripts to produce the analysis and figures for the PROPEL 16S paper: <br>
Effects of Lactobacillus reuteri supplementation on the gut microbiota in extremely preterm infants in a randomized placebo-controlled trial.<br>
https://doi.org/10.1016/j.xcrm.2021.100206

**Worklfow from raw sequences to ASV table**: <br> 
Sequences available at European Nucleotide Archive : PRJEB36531.<br>
#1.raw sequences quality check (fastqc and multiqc for raw reads) <br>
quality_check_raw.sh <br>

#2.quality trimming and filtering + quality check afterwards <br>
bbduk_trimming.sh #qtrim=l trimq=35 + qtrim=rl trimq=30 minlen=100 + quality check <br>

#3.remove phix genome NC_001422.1 from trimmed sequences + quality cehck afterwards <br>
phixgenome_removal.sh

#4.denoising and taxonomic assigment (pseudo pooling, silva database, minBoot=80)
Trimmed sequences are processed following the DADA2 Workflow with dada2 version 1.10.1 (https://benjjneb.github.io/dada2/tutorial.html)
dada2_rscript.R
dada2.sh

**Worklfow for statistical analysis**:<br>
Metadata available at https://doi.org/10.1016/j.xcrm.2021.100206 <br>
#1 Preprocessing.R <br>
#2: Figure 2: Alpha-diversity <br>
#3: Figure 3: Beta-diversity <br>
#4: Figure 4: Taxonomy <br>
#5: Figure 5: qPCR <br>
#6: Figure 6: Growth <br>


