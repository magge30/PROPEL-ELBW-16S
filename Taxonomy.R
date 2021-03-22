library(reshape2)
library(ggplot2)

# tax glom to different levels
ps3.phy<-tax_glom(ps3,"Phylum")
ps3.fam<-tax_glom(ps3,"Family")
ps3.gen<-tax_glom(ps3,"Genus")

# extract otu and tax table from phyloseq object
asv.phy = as(otu_table(ps3.phy), "matrix")
asv.fam = as(otu_table(ps3.fam), "matrix")
asv.gen = as(otu_table(ps3.gen), "matrix")

dat.phy  = as(sample_data(ps3.phy),"data.frame")
dat.fam  = as(sample_data(ps3.fam),"data.frame")
dat.gen  = as(sample_data(ps3.gen),"data.frame")

# extract taxa phyoseq object
tax.phy = as(tax_table(ps3.phy),"matrix")
tax.fam = as(tax_table(ps3.fam),"matrix")
tax.gen = as(tax_table(ps3.gen),"matrix")

tax.phy = as.data.frame(tax.phy[,"Phylum"])
tax.fam = as.data.frame(tax.fam[,"Family"])
tax.gen = as.data.frame(tax.gen[,"Genus"])

names(tax.phy)<-"Phylum"
names(tax.fam)<-"Family"
names(tax.gen)<-"Genus"

########################################################################
# Phylum
########################################################################
# arrange data for plotting
fix.data.melt<-function(asv,tax,dat){
  row.names(asv)<-tax[,1]
  df<-merge(as.data.frame(t(asv)),dat[,c("Timepoint","Supplementation")],by="row.names")
  df$nr<-paste(df$Supplementation,df$Timepoint,sep="_")
  df$nr <- gsub("Placebo_1v", 'Pl = 54', df$nr)
  df$nr <- gsub("Placebo_2v", 'Pl = 55', df$nr)
  df$nr <- gsub("Placebo_3v", 'Pl = 51', df$nr)
  df$nr <- gsub("Placebo_4v", 'Pl = 48', df$nr)
  df$nr <- gsub("Placebo_gv36", 'Pl = 41', df$nr)
  df$nr <- gsub("Placebo_2y", 'Pl = 27', df$nr)
  df$nr <- gsub("L.reuteri_1v", 'Lr = 54 ', df$nr)
  df$nr <- gsub("L.reuteri_2v", 'Lr = 54', df$nr)
  df$nr <- gsub("L.reuteri_3v", 'Lr = 51', df$nr)
  df$nr <- gsub("L.reuteri_4v", 'Lr = 53', df$nr)
  df$nr <- gsub("L.reuteri_gv36", 'Lr = 50', df$nr)
  df$nr <- gsub("L.reuteri_2y", 'Lr = 20', df$nr)
  df$SupTim<-paste(df$Supplementation,df$Timepoint,sep="_")
  dfm<-melt(df)
  #rename time point
  dfm$Timepoint <- gsub("1v","1w", dfm$Timepoint)
  dfm$Timepoint <- gsub("2v","2w", dfm$Timepoint)
  dfm$Timepoint <- gsub("3v","3w", dfm$Timepoint)
  dfm$Timepoint <- gsub("4v","4w", dfm$Timepoint)
  dfm$Timepoint <- gsub("gv36","PMW36", dfm$Timepoint)
  # to sort out facet grid by time point
  dfm$Timepoint = factor(dfm$Timepoint, levels=c("1w","2w","3w","4w","PMW36","2y"))
  # remove NA
  dfm[dfm=="NA"] <- NA 
  dfm.omit<-na.omit(dfm)
  return(dfm.omit)
}

phy.melt<-fix.data.melt(asv.phy,tax.phy,dat.phy)   # without NA

########################################################################
# Family and Genus
########################################################################
fix.data<-function(asv,tax,dat){
	row.names(asv)<-tax[,1]
	df<-merge(as.data.frame(t(asv)),dat[,c("Timepoint","Supplementation")],by="row.names")
	#rename time point
	df$Timepoint <- gsub("1v","1w", df$Timepoint)
	df$Timepoint <- gsub("2v","2w", df$Timepoint)
	df$Timepoint <- gsub("3v","3w", df$Timepoint)
	df$Timepoint <- gsub("4v","4w", df$Timepoint)
	df$Timepoint <- gsub("gv36","PMW36", df$Timepoint)
	# to sort out facet grid by time point
	df$Timepoint = factor(df$Timepoint, levels=c("1w","2w","3w","4w","PMW36","2y"))
	row.names(df)<-df$Row.names
	df<-df[,-1]
	df<-df[,colnames(df)!="NA"] 
  df<-df[, -grep("NA.*", colnames(df))]
	return(df)
}

fam<-fix.data(asv.fam,tax.fam,dat.fam)
gen<-fix.data(asv.gen,tax.gen,dat.gen)

# select only asv data
n<-ncol(fam)
fam.count<-fam[,c(1:(n-2))]
n<-ncol(gen)
gen.count<-gen[,c(1:(n-2))]

# with others as <1%
relabu1others<-function(count){
  df<-as.data.frame(t(count))
  df$relabu<-(rowSums(df)/sum(df))
  df1<-df[df$relabu>=0.01,]
  df0<-df[df$relabu<0.01,]
  df0<-as.data.frame(colSums(df0))
  names(df0)<-c("Others")
  newdf<-rbind(df1,t(df0))
  return(newdf)
}

fam.r1<-relabu1others(fam.count)
gen.r1<-relabu1others(gen.count)

# add metadata
fam1<-merge(as.data.frame(t(fam.r1)),fam[c("Timepoint","Supplementation")],by="row.names",all.y=TRUE)
gen1<-merge(as.data.frame(t(gen.r1)),gen[c("Timepoint","Supplementation")],by="row.names",all.y=TRUE)

# how much of relative abundance
gen.r2<-relabu1others(gen.count)

# Stacked bar plots
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n = 8, name = "Dark2")
brewer.pal(n = 8, name = "set2")
brewer.pal(n = 12, name = "Paired")
cbPalette<-c("#1B9E77","#FFD92F","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666","#CAB2D6")
cbPalette<-c(as.vector(brewer.pal(n = 8, name = "Dark2")), as.vector(brewer.pal(n = 12, name = "Paired")))

# Plots
p.phy<-ggplot(phy.melt, aes(fill=variable, y=value, x=Supplementation)) + 
    geom_bar(stat="identity",  position="fill") + 
    facet_grid(cols = vars(Timepoint)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5,size=5),
          axis.text.y = element_text(angle = 0, hjust = 0.5,size=8),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=10),
          legend.title=element_text(size=5),
          legend.text=element_text(size=5),
          legend.key.size=unit(0.4,"cm")) +
    labs(fill="Phylum",y="Relative abundance") +
    scale_x_discrete(limit = c("Placebo","L.reuteri"),labels=expression(Placebo, italic(L.reuteri))) +
    scale_fill_manual(values=cbPalette) 
 
p.fam<-ggplot(melt(fam1), aes(fill=variable, y=value, x=Supplementation)) + 
    geom_bar(stat="identity", position="fill") + 
    facet_grid(cols = vars(Timepoint)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5,size=5),
          axis.text.y = element_text(angle = 0, hjust = 0.5,size=8),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=10),
          legend.title=element_text(size=5),
          legend.text=element_text(size=5),
          legend.key.size=unit(0.4,"cm")) +
    labs(fill="Family",y="Relative abundance") +
    scale_x_discrete(limit = c("Placebo","L.reuteri"),labels=expression(Placebo, italic(L.reuteri))) +
    scale_fill_manual(values=cbPalette)
    
p.gen<-ggplot(melt(gen1), aes(fill=variable, y=value, x=Supplementation)) + 
    geom_bar(stat="identity", position="fill") + 
    facet_grid(cols = vars(Timepoint)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5,size=5),
          axis.text.y = element_text(angle = 0, hjust = 0.5,size=8),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=10),
          legend.title=element_text(size=5),
          legend.text=element_text(size=5),
          legend.key.size=unit(0.4,"cm")) +
    labs(fill="Genus",y="Relative abundance") +
    scale_x_discrete(limit = c("Placebo","L.reuteri"),labels=expression(Placebo, italic(L.reuteri))) +
    scale_fill_manual(values=cbPalette) 
    
ggarrange(p.phy, p.fam, p.gen, ncol = 1, nrow = 3,  labels="AUTO")

########################################################################
# Lefse
########################################################################
# Subset by timepint ps3
ps3.1w<-subset_samples(ps3,Timepoint == "1v")
ps3.2w<-subset_samples(ps3,Timepoint == "2v")
ps3.3w<-subset_samples(ps3,Timepoint == "3v")
ps3.4w<-subset_samples(ps3,Timepoint == "4v")
ps3.5w<-subset_samples(ps3,Timepoint == "gv36")

ps3.1w<- prune_taxa(taxa_sums(ps3.1w) > 0, ps3.1w) # prune taxa without observations
ps3.2w<- prune_taxa(taxa_sums(ps3.2w) > 0, ps3.2w) 
ps3.3w<- prune_taxa(taxa_sums(ps3.3w) > 0, ps3.3w) 
ps3.4w<- prune_taxa(taxa_sums(ps3.4w) > 0, ps3.4w) 
ps3.5w<- prune_taxa(taxa_sums(ps3.5w) > 0, ps3.5w) 

# tax glom to Genus,family and phylum by timepoint
ps3.1w.gen<-tax_glom(ps3.1w,"Genus")
ps3.1w.gen<-subset_taxa(ps3.1w.gen, Genus !="NA")
ps3.2w.gen<-tax_glom(ps3.2w,"Genus")
ps3.2w.gen<-subset_taxa(ps3.2w.gen, Genus !="NA")
ps3.3w.gen<-tax_glom(ps3.3w,"Genus")
ps3.3w.gen<-subset_taxa(ps3.3w.gen, Genus !="NA")
ps3.4w.gen<-tax_glom(ps3.4w,"Genus")
ps3.4w.gen<-subset_taxa(ps3.4w.gen, Genus !="NA")
ps3.5w.gen<-tax_glom(ps3.5w,"Genus")
ps3.5w.gen<-subset_taxa(ps3.5w.gen, Genus !="NA")

ps3.1w.fam<-tax_glom(ps3.1w,"Family")
ps3.1w.fam<-subset_taxa(ps3.1w.fam, Family !="NA")
ps3.2w.fam<-tax_glom(ps3.2w,"Family")
ps3.2w.fam<-subset_taxa(ps3.2w.fam, Family !="NA")
ps3.3w.fam<-tax_glom(ps3.3w,"Family")
ps3.3w.fam<-subset_taxa(ps3.3w.fam, Family !="NA")
ps3.4w.fam<-tax_glom(ps3.4w,"Family")
ps3.4w.fam<-subset_taxa(ps3.4w.fam, Family !="NA")
ps3.5w.fam<-tax_glom(ps3.5w,"Family")
ps3.5w.fam<-subset_taxa(ps3.5w.fam, Family !="NA")

ps3.1w.phy<-tax_glom(ps3.1w,"Phylum")
ps3.1w.phy<-subset_taxa(ps3.1w.phy, Phylum !="NA")

# confounders
# "Gender" 1w 3w
# "Sectio"    3w
# "ABv4.0d"      4w 
# "Location" (not confounder) 1w

# +++++ LEFSE fucntion to prepare data ++++ #
lefsedata<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# +++++ Genus Prepare data for lefse ++++ #
gen1w <- lefsedata(ps3.1w.gen)
tax_table(ps3.1w.gen)[,"Genus"]
n=nrow(gen1w)
asvname<-row.names(gen1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.1w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen1w)[3:n]<-df$noms
write.table(gen1w, "gen1w.txt",sep="\t",quote=FALSE)

gen2w <- lefsedata(ps3.2w.gen)
tax_table(ps3.2w.gen)[,"Genus"]
n=nrow(gen2w)
asvname<-row.names(gen2w)[3:n]
taxname<-as.vector(tax_table(ps3.2w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.2w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen2w)[3:n]<-df$noms
write.table(gen2w, "gen2w.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[3:n]<-df$noms
write.table(gen3w, "gen3w.txt",sep="\t",quote=FALSE)

gen4w <- lefsedata(ps3.4w.gen)
tax_table(ps3.4w.gen)[,"Genus"]
n=nrow(gen4w)
asvname<-row.names(gen4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.4w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen4w)[3:n]<-df$noms
write.table(gen4w, "gen4w.txt",sep="\t",quote=FALSE)

gen5w <- lefsedata(ps3.5w.gen)
tax_table(ps3.5w.gen)[,"Genus"]
n=nrow(gen5w)
asvname<-row.names(gen5w)[3:n]
taxname<-as.vector(tax_table(ps3.5w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.5w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen5w)[3:n]<-df$noms
write.table(gen5w, "gen5w.txt",sep="\t",quote=FALSE)

# +++++ Family Prepare data for lefse ++++ #
fam1w <- lefsedata(ps3.1w.fam)
tax_table(ps3.1w.fam)[,"Family"]
n=nrow(fam1w)
asvname<-row.names(fam1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.1w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam1w)[3:n]<-df$noms
write.table(fam1w, "fam1w.txt",sep="\t",quote=FALSE)

fam2w <- lefsedata(ps3.2w.fam)
tax_table(ps3.2w.fam)[,"Family"]
n=nrow(fam2w)
asvname<-row.names(fam2w)[3:n]
taxname<-as.vector(tax_table(ps3.2w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.2w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam2w)[3:n]<-df$noms
write.table(fam2w, "fam2w.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[3:n]<-df$noms
write.table(fam3w, "fam3w.txt",sep="\t",quote=FALSE)

fam4w <- lefsedata(ps3.4w.fam)
tax_table(ps3.4w.fam)[,"Family"]
n=nrow(fam4w)
asvname<-row.names(fam4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.4w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam4w)[3:n]<-df$noms
write.table(fam4w, "fam4w.txt",sep="\t",quote=FALSE)

fam5w <- lefsedata(ps3.5w.fam)
tax_table(ps3.5w.fam)[,"Family"]
n=nrow(fam5w)
asvname<-row.names(fam5w)[3:n]
taxname<-as.vector(tax_table(ps3.5w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.5w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam5w)[3:n]<-df$noms
write.table(fam5w, "fam5w.txt",sep="\t",quote=FALSE)

# +++++ Phylum Prepare data for lefse ++++ #
phy1w <- lefsedata(ps3.1w.phy)
tax_table(ps3.1w.phy)[,"Phylum"]
n=nrow(phy1w)
asvname<-row.names(phy1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.phy)[,"Phylum"])
taxasv<-row.names(tax_table(ps3.1w.phy))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(phy1w)[3:n]<-df$noms
write.table(phy1w, "phy1w.txt",sep="\t",quote=FALSE)


##########################################################
# GALAXY HUTLAB
##########################################################

# GO TO galaxy: https://huttenhower.sph.harvard.edu/galaxy/

#1 - Get Data: define Type "Tabular" (/Users/magge30/Documents/Research/Postdoc_2017_ThomasAbrahamsson/PROPEL/16S_DataAnalysis/R_analyis/R_results/LEFSE/lefse_manuscript)
#2 - LEfSe 

# Input_MS_DSM17938: includes DSM and lactobacillus
# workflow: MS_0.01
# A) Format data for lefe : Class: "Supplementation", Subject: "SampleID", Subclass: confounders
# B) LDA Effect Size: alpha value: 0.01, LDS score 2
# C) Plot leefse results
# D) 

########################################################################
# CONFOUNDERS
########################################################################

#++++ CONFOUNDERS ++++#

# 1w confounder gender
lefsedata1<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation","Gender")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 1w only Gender
lefsedata1g<-function(ps){ 
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Gender")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 3w confounder gender
lefsedata3g<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation","Gender")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 3w confounder sectio
lefsedata3d<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation","Sectio")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 3w  only sectio
lefsedata3s<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Sectio")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 3w  only gender
lefsedata3gc<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Gender")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 4w confounderr antibiotic
lefsedata4<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation","ABv4.0d")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 4w only antibiotic
lefsedata4ac<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","ABv4.0d")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

phy1w <- lefsedata(ps3.1w.phy)
tax_table(ps3.1w.phy)[,"Phylum"]
n=nrow(phy1w)
asvname<-row.names(phy1w)[4:n]
taxname<-as.vector(tax_table(ps3.1w.phy)[,"Phylum"])
taxasv<-row.names(tax_table(ps3.1w.phy))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(phy1w)[4:n]<-df$noms
write.table(phy1w, "phy1w.txt",sep="\t",quote=FALSE)

phy1w <- lefsedata1g(ps3.1w.phy)
tax_table(ps3.1w.phy)[,"Phylum"]
n=nrow(phy1w)
asvname<-row.names(phy1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.phy)[,"Phylum"])
taxasv<-row.names(tax_table(ps3.1w.phy))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(phy1w)[3:n]<-df$noms
write.table(phy1w, "phy1wg.txt",sep="\t",quote=FALSE)

fam1w <- lefsedata1(ps3.1w.fam)
tax_table(ps3.1w.fam)[,"Family"]
n=nrow(fam1w)
asvname<-row.names(fam1w)[4:n]
taxname<-as.vector(tax_table(ps3.1w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.1w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam1w)[4:n]<-df$noms
write.table(fam1w, "fam1w.txt",sep="\t",quote=FALSE)

fam1w <- lefsedata1g(ps3.1w.fam)
tax_table(ps3.1w.fam)[,"Family"]
n=nrow(fam1w)
asvname<-row.names(fam1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.1w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam1w)[3:n]<-df$noms
write.table(fam1w, "fam1wg.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata3g(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[4:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[4:n]<-df$noms
write.table(fam3w, "fam3wg.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata3d(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[4:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[4:n]<-df$noms
write.table(fam3w, "fam3wd.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata3s(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[3:n]<-df$noms
write.table(fam3w, "fam3ws.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata3gc(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[3:n]<-df$noms
write.table(fam3w, "fam3wgc.txt",sep="\t",quote=FALSE)

fam4w <- lefsedata4(ps3.4w.fam)
tax_table(ps3.4w.fam)[,"Family"]
n=nrow(fam4w)
asvname<-row.names(fam4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.4w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam4w)[3:n]<-df$noms
write.table(fam4w, "fam4w.txt",sep="\t",quote=FALSE)

fam4w <- lefsedata4ac(ps3.4w.fam)
tax_table(ps3.4w.fam)[,"Family"]
n=nrow(fam4w)
asvname<-row.names(fam4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.4w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam4w)[3:n]<-df$noms
write.table(fam4w, "fam4wac.txt",sep="\t",quote=FALSE)

gen1w <- lefsedata1(ps3.1w.gen)
tax_table(ps3.1w.gen)[,"Genus"]
n=nrow(gen1w)
asvname<-row.names(gen1w)[4:n]
taxname<-as.vector(tax_table(ps3.1w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.1w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen1w)[4:n]<-df$noms
write.table(gen1w, "gen1w.txt",sep="\t",quote=FALSE)

gen1w <- lefsedata1g(ps3.1w.gen)
tax_table(ps3.1w.gen)[,"Genus"]
n=nrow(gen1w)
asvname<-row.names(gen1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.1w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen1w)[3:n]<-df$noms
write.table(gen1w, "gen1wg.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata3g(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[4:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[4:n]<-df$noms
write.table(gen3w, "gen3wg.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata3d(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[4:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[4:n]<-df$noms
write.table(gen3w, "gen3wd.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata3s(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[3:n]<-df$noms
write.table(gen3w, "gen3ws.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata3gc(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[3:n]<-df$noms
write.table(gen3w, "gen3wgc.txt",sep="\t",quote=FALSE)

gen4w <- lefsedata4(ps3.4w.gen)
tax_table(ps3.4w.gen)[,"Genus"]
n=nrow(gen4w)
asvname<-row.names(gen4w)[4:n]
taxname<-as.vector(tax_table(ps3.4w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.4w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen4w)[4:n]<-df$noms
write.table(gen4w, "gen4w.txt",sep="\t",quote=FALSE)

gen4w <- lefsedata4ac(ps3.4w.gen)
tax_table(ps3.4w.gen)[,"Genus"]
n=nrow(gen4w)
asvname<-row.names(gen4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.4w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen4w)[3:n]<-df$noms
write.table(gen4w, "gen4wac.txt",sep="\t",quote=FALSE)

########################################################################
# Lefse
########################################################################

# tax glom to Genus,family and phylum by timepoint
ps3.1w.gen<-tax_glom(ps3.1w,"Genus")
ps3.1w.gen<-subset_taxa(ps3.1w.gen, Genus !="NA")
ps3.2w.gen<-tax_glom(ps3.2w,"Genus")
ps3.2w.gen<-subset_taxa(ps3.2w.gen, Genus !="NA")
ps3.3w.gen<-tax_glom(ps3.3w,"Genus")
ps3.3w.gen<-subset_taxa(ps3.3w.gen, Genus !="NA")
ps3.4w.gen<-tax_glom(ps3.4w,"Genus")
ps3.4w.gen<-subset_taxa(ps3.4w.gen, Genus !="NA")
ps3.5w.gen<-tax_glom(ps3.5w,"Genus")
ps3.5w.gen<-subset_taxa(ps3.5w.gen, Genus !="NA")

ps3.1w.fam<-tax_glom(ps3.1w,"Family")
ps3.1w.fam<-subset_taxa(ps3.1w.fam, Family !="NA")
ps3.2w.fam<-tax_glom(ps3.2w,"Family")
ps3.2w.fam<-subset_taxa(ps3.2w.fam, Family !="NA")
ps3.3w.fam<-tax_glom(ps3.3w,"Family")
ps3.3w.fam<-subset_taxa(ps3.3w.fam, Family !="NA")
ps3.4w.fam<-tax_glom(ps3.4w,"Family")
ps3.4w.fam<-subset_taxa(ps3.4w.fam, Family !="NA")
ps3.5w.fam<-tax_glom(ps3.5w,"Family")
ps3.5w.fam<-subset_taxa(ps3.5w.fam, Family !="NA")

ps3.1w.phy<-tax_glom(ps3.1w,"Phylum")
ps3.1w.phy<-subset_taxa(ps3.1w.phy, Phylum !="NA")

# confounders
# "Gender" 1w 3w
# "Sectio"    3w
# "ABv4.0d"      4w 
# "Location" (not confounder) 1w

# +++++ LEFSE fucntion to prepare data ++++ #
lefsedata<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# +++++ Genus Prepare data for lefse ++++ #
gen1w <- lefsedata(ps3.1w.gen)
tax_table(ps3.1w.gen)[,"Genus"]
n=nrow(gen1w)
asvname<-row.names(gen1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.1w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen1w)[3:n]<-df$noms
write.table(gen1w, "gen1w.txt",sep="\t",quote=FALSE)

gen2w <- lefsedata(ps3.2w.gen)
tax_table(ps3.2w.gen)[,"Genus"]
n=nrow(gen2w)
asvname<-row.names(gen2w)[3:n]
taxname<-as.vector(tax_table(ps3.2w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.2w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen2w)[3:n]<-df$noms
write.table(gen2w, "gen2w.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[3:n]<-df$noms
write.table(gen3w, "gen3w.txt",sep="\t",quote=FALSE)

gen4w <- lefsedata(ps3.4w.gen)
tax_table(ps3.4w.gen)[,"Genus"]
n=nrow(gen4w)
asvname<-row.names(gen4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.4w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen4w)[3:n]<-df$noms
write.table(gen4w, "gen4w.txt",sep="\t",quote=FALSE)

gen5w <- lefsedata(ps3.5w.gen)
tax_table(ps3.5w.gen)[,"Genus"]
n=nrow(gen5w)
asvname<-row.names(gen5w)[3:n]
taxname<-as.vector(tax_table(ps3.5w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.5w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen5w)[3:n]<-df$noms
write.table(gen5w, "gen5w.txt",sep="\t",quote=FALSE)

# +++++ Family Prepare data for lefse ++++ #
fam1w <- lefsedata(ps3.1w.fam)
tax_table(ps3.1w.fam)[,"Family"]
n=nrow(fam1w)
asvname<-row.names(fam1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.1w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam1w)[3:n]<-df$noms
write.table(fam1w, "fam1w.txt",sep="\t",quote=FALSE)

fam2w <- lefsedata(ps3.2w.fam)
tax_table(ps3.2w.fam)[,"Family"]
n=nrow(fam2w)
asvname<-row.names(fam2w)[3:n]
taxname<-as.vector(tax_table(ps3.2w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.2w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam2w)[3:n]<-df$noms
write.table(fam2w, "fam2w.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[3:n]<-df$noms
write.table(fam3w, "fam3w.txt",sep="\t",quote=FALSE)

fam4w <- lefsedata(ps3.4w.fam)
tax_table(ps3.4w.fam)[,"Family"]
n=nrow(fam4w)
asvname<-row.names(fam4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.4w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam4w)[3:n]<-df$noms
write.table(fam4w, "fam4w.txt",sep="\t",quote=FALSE)

fam5w <- lefsedata(ps3.5w.fam)
tax_table(ps3.5w.fam)[,"Family"]
n=nrow(fam5w)
asvname<-row.names(fam5w)[3:n]
taxname<-as.vector(tax_table(ps3.5w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.5w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam5w)[3:n]<-df$noms
write.table(fam5w, "fam5w.txt",sep="\t",quote=FALSE)

# +++++ Phylum Prepare data for lefse ++++ #
phy1w <- lefsedata(ps3.1w.phy)
tax_table(ps3.1w.phy)[,"Phylum"]
n=nrow(phy1w)
asvname<-row.names(phy1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.phy)[,"Phylum"])
taxasv<-row.names(tax_table(ps3.1w.phy))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(phy1w)[3:n]<-df$noms
write.table(phy1w, "phy1w.txt",sep="\t",quote=FALSE)


# GO TO galaxy: https://huttenhower.sph.harvard.edu/galaxy/
#1 - Get Data: define Type "Tabular" 
#2 - LEfSe 

# workflow: MS_0.01
# A) Format data for lefe : Class: "Supplementation", Subject: "SampleID", Subclass: confounders
# B) LDA Effect Size: alpha value: 0.01, LDS score 2
# C) Plot leefse results
# D) 

#++++ CONFOUNDERS ++++#

# 1w confounder gender
lefsedata1<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation","Gender")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 1w only Gender
lefsedata1g<-function(ps){ 
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Gender")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 3w confounder gender
lefsedata3g<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation","Gender")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 3w confounder sectio
lefsedata3d<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation","Sectio")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 3w  only sectio
lefsedata3s<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Sectio")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 3w  only gender
lefsedata3gc<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Gender")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 4w confounderr antibiotic
lefsedata4<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","Supplementation","ABv4.0d")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

# 4w only antibiotic
lefsedata4ac<-function(ps){
  asv = as(otu_table(ps), "matrix")
  dat = as(sample_data(ps),"data.frame")
  dat = dat[,c("SampleID","ABv4.0d")]
  df<-merge(dat,t(asv),by="row.names")
  row.names(df)<-df$Row.names
  df<-as.data.frame(t(df))
  df<-df[-1,]
  names(df)<-NULL
  return(df) 
}

phy1w <- lefsedata(ps3.1w.phy)
tax_table(ps3.1w.phy)[,"Phylum"]
n=nrow(phy1w)
asvname<-row.names(phy1w)[4:n]
taxname<-as.vector(tax_table(ps3.1w.phy)[,"Phylum"])
taxasv<-row.names(tax_table(ps3.1w.phy))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(phy1w)[4:n]<-df$noms
write.table(phy1w, "phy1w.txt",sep="\t",quote=FALSE)

phy1w <- lefsedata1g(ps3.1w.phy)
tax_table(ps3.1w.phy)[,"Phylum"]
n=nrow(phy1w)
asvname<-row.names(phy1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.phy)[,"Phylum"])
taxasv<-row.names(tax_table(ps3.1w.phy))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(phy1w)[3:n]<-df$noms
write.table(phy1w, "phy1wg.txt",sep="\t",quote=FALSE)

fam1w <- lefsedata1(ps3.1w.fam)
tax_table(ps3.1w.fam)[,"Family"]
n=nrow(fam1w)
asvname<-row.names(fam1w)[4:n]
taxname<-as.vector(tax_table(ps3.1w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.1w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam1w)[4:n]<-df$noms
write.table(fam1w, "fam1w.txt",sep="\t",quote=FALSE)

fam1w <- lefsedata1g(ps3.1w.fam)
tax_table(ps3.1w.fam)[,"Family"]
n=nrow(fam1w)
asvname<-row.names(fam1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.1w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam1w)[3:n]<-df$noms
write.table(fam1w, "fam1wg.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata3g(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[4:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[4:n]<-df$noms
write.table(fam3w, "fam3wg.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata3d(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[4:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[4:n]<-df$noms
write.table(fam3w, "fam3wd.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata3s(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[3:n]<-df$noms
write.table(fam3w, "fam3ws.txt",sep="\t",quote=FALSE)

fam3w <- lefsedata3gc(ps3.3w.fam)
tax_table(ps3.3w.fam)[,"Family"]
n=nrow(fam3w)
asvname<-row.names(fam3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.3w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam3w)[3:n]<-df$noms
write.table(fam3w, "fam3wgc.txt",sep="\t",quote=FALSE)

fam4w <- lefsedata4(ps3.4w.fam)
tax_table(ps3.4w.fam)[,"Family"]
n=nrow(fam4w)
asvname<-row.names(fam4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.4w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam4w)[3:n]<-df$noms
write.table(fam4w, "fam4w.txt",sep="\t",quote=FALSE)

fam4w <- lefsedata4ac(ps3.4w.fam)
tax_table(ps3.4w.fam)[,"Family"]
n=nrow(fam4w)
asvname<-row.names(fam4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.fam)[,"Family"])
taxasv<-row.names(tax_table(ps3.4w.fam))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(fam4w)[3:n]<-df$noms
write.table(fam4w, "fam4wac.txt",sep="\t",quote=FALSE)

gen1w <- lefsedata1(ps3.1w.gen)
tax_table(ps3.1w.gen)[,"Genus"]
n=nrow(gen1w)
asvname<-row.names(gen1w)[4:n]
taxname<-as.vector(tax_table(ps3.1w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.1w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen1w)[4:n]<-df$noms
write.table(gen1w, "gen1w.txt",sep="\t",quote=FALSE)

gen1w <- lefsedata1g(ps3.1w.gen)
tax_table(ps3.1w.gen)[,"Genus"]
n=nrow(gen1w)
asvname<-row.names(gen1w)[3:n]
taxname<-as.vector(tax_table(ps3.1w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.1w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen1w)[3:n]<-df$noms
write.table(gen1w, "gen1wg.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata3g(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[4:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[4:n]<-df$noms
write.table(gen3w, "gen3wg.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata3d(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[4:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[4:n]<-df$noms
write.table(gen3w, "gen3wd.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata3s(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[3:n]<-df$noms
write.table(gen3w, "gen3ws.txt",sep="\t",quote=FALSE)

gen3w <- lefsedata3gc(ps3.3w.gen)
tax_table(ps3.3w.gen)[,"Genus"]
n=nrow(gen3w)
asvname<-row.names(gen3w)[3:n]
taxname<-as.vector(tax_table(ps3.3w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.3w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen3w)[3:n]<-df$noms
write.table(gen3w, "gen3wgc.txt",sep="\t",quote=FALSE)

gen4w <- lefsedata4(ps3.4w.gen)
tax_table(ps3.4w.gen)[,"Genus"]
n=nrow(gen4w)
asvname<-row.names(gen4w)[4:n]
taxname<-as.vector(tax_table(ps3.4w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.4w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen4w)[4:n]<-df$noms
write.table(gen4w, "gen4w.txt",sep="\t",quote=FALSE)

gen4w <- lefsedata4ac(ps3.4w.gen)
tax_table(ps3.4w.gen)[,"Genus"]
n=nrow(gen4w)
asvname<-row.names(gen4w)[3:n]
taxname<-as.vector(tax_table(ps3.4w.gen)[,"Genus"])
taxasv<-row.names(tax_table(ps3.4w.gen))
df<-data.frame(asvname,taxasv,taxname)
df$noms<-paste(df$taxname,"_",df$asvname)
row.names(gen4w)[3:n]<-df$noms
write.table(gen4w, "gen4wac.txt",sep="\t",quote=FALSE)

