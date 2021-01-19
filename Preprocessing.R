# packages
library("phyloseq")
library("vegan")
library("plyr")

# load data: metadata, count_tab and tax_tab

# phyloseq object 
ps <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
               sample_data(data), 
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