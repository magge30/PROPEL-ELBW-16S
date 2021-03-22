library(dada2); packageVersion("dada2")

# set working directory
setwd("working/directory")
#list.files()

#path to the directory containing the fastq files
path <- "working/directory" 
#list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_trimmed_phixunmatch.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_trimmed_phixunmatch.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
 
# Assign the filenames for the filtered fastq.gz files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

# remove quality trim/filter: default settings
filtered_out<- filterAndTrim(fnFs, filtFs, fnRs, filtRs)                                                      
dim(filtered_out)
head(filtered_out)

# learn error rates 
errF <- learnErrors(filtFs, multithread=FALSE) 
errR <- learnErrors(filtRs, multithread=FALSE)

# dereplication 
derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)
                                                         
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

######################################################################################################
# infer ASVs on both forward and reverse reads by pseudo pooling
######################################################################################################
dadaFs.pseudo <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs.pseudo <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")
 
#merge paired
mergers.pseudo <- mergePairs(dadaFs.pseudo, derepFs, dadaRs.pseudo, derepRs, trimOverhang=TRUE, verbose=TRUE)

# construct a amplicon sequence variant table (ASV) table
seqtab.pseudo <- makeSequenceTable(mergers.pseudo)
dim(seqtab.pseudo)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.pseudo)))
write.table(table(nchar(getSequences(seqtab.pseudo))),"seqlength_pseudo.txt",sep="\t") 

# To remove non-targe-length sequences
seqtab2.pseudo <- seqtab.pseudo[,nchar(colnames(seqtab.pseudo)) %in% seq(100,500)]
seqtab.pseudo<-seqtab2.pseudo
write.table(table(nchar(getSequences(seqtab2.pseudo))),"seqlenght2_pseudo.txt",sep="\t")

# chimeras removal
seqtab.nochim.pseudo <- removeBimeraDenovo(seqtab.pseudo, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim.pseudo)

# how much in terms of abundance did the chimera hold:
sum(seqtab.nochim.pseudo)/sum(seqtab.pseudo)

# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
summary_tab.pseudo <- data.frame(row.names=sample.names, dada2_input=filtered_out[,1], filtered=filtered_out[,2], dada_f=sapply(dadaFs.pseudo, getN), dada_r=sapply(dadaRs.pseudo, getN), merged=sapply(mergers.pseudo, getN), nonchim=rowSums(seqtab.nochim.pseudo), final_perc_reads_retained=round(rowSums(seqtab.nochim.pseudo)/filtered_out[,1]*100, 1))
write.table(summary_tab.pseudo,"summary_tab_pseudo.txt",sep="\t")

# ASSING TAXONOMY
#set.seed(100) # Initialize random number generator for reproducibility
taxa50.pseudo <- assignTaxonomy(seqtab.nochim.pseudo, "/working/directory/silva_nr_v132_train_set.fa.gz", multithread=T, tryRC=T,minBoot=50)
taxa80.pseudo <- assignTaxonomy(seqtab.nochim.pseudo, "/working/directory/silva_nr_v132_train_set.fa.gz", multithread=T, tryRC=T,minBoot=80)

# add species 
taxa50.pseudo.species <- addSpecies(taxa50.pseudo, "/working/directory/silva_species_assignment_v132.fa.gz")
#unname(taxa50.pseudo.species)

taxa50.pseudo.speciesM <- addSpecies(taxa50.pseudo, "/working/directory/silva_species_assignment_v132.fa.gz",allowMultiple=TRUE)
#unname(taxa50.pseudo.speciesM)

taxa80.pseudo.species <- addSpecies(taxa80.pseudo, "/working/directory/silva_species_assignment_v132.fa.gz")
#unname(taxa80.pseudo.species)

taxa80.pseudo.speciesM <- addSpecies(taxa80.pseudo, "/working/directory/silva/silva_species_assignment_v132.fa.gz",allowMultiple=TRUE)
#unname(taxa80.pseudo.speciesM)
                                                        
# EXTRACTING STANDARD OUTPUTS: fasta file, count table and taxonmy table
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs.pseudo <- colnames(seqtab.nochim.pseudo)
asv_headers.pseudo <- vector(dim(seqtab.nochim.pseudo)[2], mode="character")

for (i in 1:dim(seqtab.nochim.pseudo)[2]) {
  asv_headers.pseudo[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta.pseudo <- c(rbind(asv_headers.pseudo, asv_seqs.pseudo))
write(asv_fasta.pseudo, "ASVs_pseudo.fa")

# count table:
asv_tab.pseudo <- t(seqtab.nochim.pseudo)
row.names(asv_tab.pseudo) <- sub(">", "", asv_headers.pseudo)
write.table(asv_tab.pseudo, "ASVs_counts_pseudo.txt", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax50.pseudo <- taxa50.pseudo
row.names(asv_tax50.pseudo) <- sub(">", "", asv_headers.pseudo)
write.table(asv_tax50.pseudo, "ASVs_taxonomy50_pseudo.txt", sep="\t", quote=F, col.names=NA)

asv_tax80.pseudo <- taxa80.pseudo
row.names(asv_tax80.pseudo) <- sub(">", "", asv_headers.pseudo)
write.table(asv_tax80.pseudo, "ASVs_taxonomy80_pseudo.txt", sep="\t", quote=F, col.names=NA)

asv_tax50.pseudo.speciesM <- taxa50.pseudo.speciesM
row.names(asv_tax50.pseudo.speciesM) <- sub(">", "", asv_headers.pseudo)
write.table(asv_tax50.pseudo.speciesM, "ASVs_taxonomy50_pseudo_spM.txt", sep="\t", quote=F, col.names=NA)

asv_tax80.pseudo.speciesM <- taxa80.pseudo.speciesM
row.names(asv_tax80.pseudo.speciesM) <- sub(">", "", asv_headers.pseudo)
write.table(asv_tax80.pseudo.speciesM, "ASVs_taxonomy80_pseudo_spM.txt", sep="\t", quote=F, col.names=NA)
