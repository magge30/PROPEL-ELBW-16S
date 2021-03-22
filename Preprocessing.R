# packages
library(phyloseq)
library(vegan)
library(plyr)

# required input data: metadata (Table S10), count_tab and tax_tab (output of the "Worklfow from raw sequences to ASV table")

# phyloseq object 
ps <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
               sample_data(metadata), 
               tax_table(tax_tab))

# Remove mock sample from phyloseq
ps <- prune_samples(sample_names(ps) != "mock1", ps) 
ps <- prune_samples(sample_names(ps) != "mock2", ps)
ps <- prune_samples(sample_names(ps) != "mock3", ps)
ps <- prune_samples(sample_names(ps) != "mock4", ps)
ps <- prune_samples(sample_names(ps) != "mock5", ps)
ps <- prune_samples(sample_names(ps) != "mock6", ps)
ps <- prune_samples(sample_names(ps) != "mock7", ps)
ps <- prune_samples(sample_names(ps) != "mock8", ps)
ps <- prune_samples(sample_names(ps) != "mock9", ps)

# prune taxa without observations after removing the mock samples
ps<- prune_taxa(taxa_sums(ps) > 0, ps)

## preprocessing
# ps1 - Taxa
ps1<-subset_taxa(ps, Kingdom == "Bacteria" & 
                    Family!= "mitochondria" &
                    Phylum != "Cyanobacteria") # Order != "Chloroplast" &

# ps2 - Prevalence
prev1 = apply(X = otu_table(ps1),
                 MARGIN = 1,
                 FUN = function(x){sum(x > 0)})

prev1 = data.frame(Prevalence = prev1,
                      TotalAbundance = taxa_sums(ps1),
                      tax_table(ps1))
keepTaxa = rownames(prev1)[(prev1$Prevalence >= 1) & (prev1$TotalAbundance > 30)]
ps2 = prune_taxa(keepTaxa, ps1)

# ps3 - Rarefy (241-1v, L.reuteri and 243-3v Placebo)
asv = as(otu_table(ps2), "matrix")
asvt<-t(asv)
set.seed(3683)
asv.rare<-rrarefy(asvt,134932)
sort(rowSums(asv.rare))

ps3 <- phyloseq(otu_table(t(asv.rare), taxa_are_rows=TRUE), 
               sample_data(data), 
               tax_table(tax_tab))

ps3<- prune_taxa(taxa_sums(ps3) > 0, ps3)

# Subset by timepoint
ps3.1w<-subset_samples(ps3,Timepoint == "1v")
ps3.2w<-subset_samples(ps3,Timepoint == "2v")
ps3.3w<-subset_samples(ps3,Timepoint == "3v")
ps3.4w<-subset_samples(ps3,Timepoint == "4v")
ps3.gw36<-subset_samples(ps3,Timepoint == "gv36")
ps3.2y<-subset_samples(ps3,Timepoint == "2y")

ps3.1w<- prune_taxa(taxa_sums(ps3.1w) > 0, ps3.1w) # prune taxa without observations
ps3.2w<- prune_taxa(taxa_sums(ps3.2w) > 0, ps3.2w) 
ps3.3w<- prune_taxa(taxa_sums(ps3.3w) > 0, ps3.3w) 
ps3.4w<- prune_taxa(taxa_sums(ps3.4w) > 0, ps3.4w) 
ps3.gw36<- prune_taxa(taxa_sums(ps3.gw36) > 0, ps3.gw36) 
ps3.2y<- prune_taxa(taxa_sums(ps3.2y) > 0, ps3.2y)