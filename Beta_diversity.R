library(DESeq2)
library(ggplot2)
library(ggpubr)

## normalization ##
#1w
dds = phyloseq_to_deseq2(ps3.1w, ~Supplementation)
deseq_counts <- estimateSizeFactors(dds, type = "poscounts") # add pseudcounts
NormVST = t(assay(varianceStabilizingTransformation(deseq_counts))) # normalization vst
ps3.1w.VST = ps3.1w
otu_table(ps3.1w.VST)<-otu_table(NormVST,FALSE)

#2w
dds = phyloseq_to_deseq2(ps3.2w, ~Supplementation)
deseq_counts <- estimateSizeFactors(dds, type = "poscounts") # add pseudcounts
NormVST = t(assay(varianceStabilizingTransformation(deseq_counts))) # normalization vst
ps3.2w.VST = ps3.2w
otu_table(ps3.2w.VST)<-otu_table(NormVST,FALSE)

#3w
dds = phyloseq_to_deseq2(ps3.3w, ~Supplementation)
deseq_counts <- estimateSizeFactors(dds, type = "poscounts") # add pseudcounts
NormVST = t(assay(varianceStabilizingTransformation(deseq_counts))) # normalization vst
ps3.3w.VST = ps3.3w
otu_table(ps3.3w.VST)<-otu_table(NormVST,FALSE)

#4w
dds = phyloseq_to_deseq2(ps3.4w, ~Supplementation)
deseq_counts <- estimateSizeFactors(dds, type = "poscounts") # add pseudcounts
NormVST = t(assay(varianceStabilizingTransformation(deseq_counts))) # normalization vst
ps3.4w.VST = ps3.4w
otu_table(ps3.4w.VST)<-otu_table(NormVST,FALSE)

#gw36
dds = phyloseq_to_deseq2(ps3.gw36, ~Supplementation)
deseq_counts <- estimateSizeFactors(dds, type = "poscounts") # add pseudcounts
NormVST = t(assay(varianceStabilizingTransformation(deseq_counts))) # normalization vst
ps3.gw36.VST = ps3.gw36
otu_table(ps3.gw36.VST)<-otu_table(NormVST,FALSE)

#2y
dds = phyloseq_to_deseq2(ps3.2y, ~Supplementation)
deseq_counts <- estimateSizeFactors(dds, type = "poscounts") # add pseudcounts
NormVST = t(assay(varianceStabilizingTransformation(deseq_counts))) # normalization vst
ps3.2y.VST = ps3.2y
otu_table(ps3.2y.VST)<-otu_table(NormVST,FALSE)

## nmds + envfit ##
#1w
otu<-as(otu_table(ps3.1w.VST),"matrix")
otu<-as.data.frame(otu)

tax=as(tax_table(ps3.1w.VST),"matrix")
tax<-as.data.frame(tax)

dat<-as(sample_data(ps3.1w.VST),"data.frame")

nmds.1w <- ordinate(ps3.1w.VST, method="NMDS", distance="euclidean") 
scrs <- as.data.frame(scores(nmds.1w, display = "sites"))
scrs <- cbind(scrs,Supplementation = dat$Supplementation, Location=dat$Location, Antibiotic=dat$Antibiotic)
scrs$Location<-as.character(scrs$Location)
scrs$Location<-gsub("0","Linköping",scrs$Location) # 0 is L and Location is Inclusion_site in Table S10
scrs$Location<-gsub("1","Stockholm",scrs$Location) # 1 is S and Location is Inclusion_site in Table S10 
scrs$locAnt<-paste(scrs$Location,scrs$Antibiotic,sep="_")
set.seed(3687)
spec <- envfit(nmds.1w, otu, perm=999)

spp.scrs <- as.data.frame(scores(spec, display="vectors"))
spp.scrs <- cbind(spp.scrs, Species =rownames(spp.scrs))
spp.scrs <- cbind(spp.scrs, pVal=spec$vectors$pvals)
spp.scrs <- cbind(spp.scrs, r=spec$vectors[2])

spp.scrs <- cbind(spp.scrs, Taxa=tax$Genus)
spp.scrs<-subset(spp.scrs,Taxa!="NA")
spp.sig<-subset(spp.scrs,pVal<=0.01)
spp.sig<-subset(spp.sig,r>=0.3)
spp.sig[,c(1,2)]<-spp.sig[,c(1,2)]*50
spp.plot<-spp.sig[c("ASV_115","ASV_8","ASV_1"),]

scrs1=scrs
spp1=spp.plot

#2w
otu<-as(otu_table(ps3.2w.VST),"matrix")
otu<-as.data.frame(otu)

tax=as(tax_table(ps3.2w.VST),"matrix")
tax<-as.data.frame(tax)

dat<-as(sample_data(ps3.2w.VST),"data.frame")

nmds.2w <- ordinate(ps3.2w.VST, method="NMDS", distance="euclidean") 
scrs <- as.data.frame(scores(nmds.2w, display = "sites"))
scrs <- cbind(scrs,Supplementation = dat$Supplementation, Location=dat$Location)
scrs$Location<-as.character(scrs$Location)
scrs$Location<-gsub("0","Linköping",scrs$Location)
scrs$Location<-gsub("1","Stockholm",scrs$Location)
set.seed(3687)
spec <- envfit(nmds.2w, otu, perm=999)

spp.scrs <- as.data.frame(scores(spec, display="vectors"))
spp.scrs <- cbind(spp.scrs, Species =rownames(spp.scrs))
spp.scrs <- cbind(spp.scrs, pVal=spec$vectors$pvals)
spp.scrs <- cbind(spp.scrs, r=spec$vectors[2])

spp.scrs <- cbind(spp.scrs, Taxa=tax$Genus)
spp.scrs<-subset(spp.scrs,Taxa!="NA")
spp.sig<-subset(spp.scrs,pVal<=0.01)
spp.sig<-subset(spp.sig,r>=0.3)
spp.sig[,c(1,2)]<-spp.sig[,c(1,2)]*0.02
spp.plot<-spp.sig[c("ASV_17","ASV_469"),]

scrs2=scrs
spp2=spp.plot

#3w
otu<-as(otu_table(ps3.3w.VST),"matrix")
otu<-as.data.frame(otu)

tax=as(tax_table(ps3.3w.VST),"matrix")
tax<-as.data.frame(tax)

dat<-as(sample_data(ps3.3w.VST),"data.frame")

nmds.3w <- ordinate(ps3.3w.VST, method="NMDS", distance="euclidean") 
scrs <- as.data.frame(scores(nmds.3w, display = "sites"))
scrs <- cbind(scrs,Supplementation = dat$Supplementation, Location=dat$Location)
scrs$Location<-as.character(scrs$Location)
scrs$Location<-gsub("0","Linköping",scrs$Location)
scrs$Location<-gsub("1","Stockholm",scrs$Location)
set.seed(3687)
spec <- envfit(nmds.3w, otu, perm=999)

spp.scrs <- as.data.frame(scores(spec, display="vectors"))
spp.scrs <- cbind(spp.scrs, Species =rownames(spp.scrs))
spp.scrs <- cbind(spp.scrs, pVal=spec$vectors$pvals)
spp.scrs <- cbind(spp.scrs, r=spec$vectors[2])

spp.scrs <- cbind(spp.scrs, Taxa=tax$Genus)
spp.scrs<-subset(spp.scrs,Taxa!="NA")
spp.sig<-subset(spp.scrs,pVal<=0.01)
spp.sig<-subset(spp.sig,r>=0.3)
spp.sig[,c(1,2)]<-spp.sig[,c(1,2)]*50
spp.plot<-spp.sig[c("ASV_15","ASV_5","ASV_27"),]

scrs3=scrs
spp3=spp.plot

#4w
otu<-as(otu_table(ps3.4w.VST),"matrix")
otu<-as.data.frame(otu)

tax=as(tax_table(ps3.4w.VST),"matrix")
tax<-as.data.frame(tax)

dat<-as(sample_data(ps3.4w.VST),"data.frame")

nmds.4w <- ordinate(ps3.4w.VST, method="NMDS", distance="euclidean") 
scrs <- as.data.frame(scores(nmds.4w, display = "sites"))
scrs <- cbind(scrs,Supplementation = dat$Supplementation, Location=dat$Location)
scrs$Location<-as.character(scrs$Location)
scrs$Location<-gsub("0","Linköping",scrs$Location)
scrs$Location<-gsub("1","Stockholm",scrs$Location)
set.seed(3687)
spec <- envfit(nmds.4w, otu, perm=999)

spp.scrs <- as.data.frame(scores(spec, display="vectors"))
spp.scrs <- cbind(spp.scrs, Species =rownames(spp.scrs))
spp.scrs <- cbind(spp.scrs, pVal=spec$vectors$pvals)
spp.scrs <- cbind(spp.scrs, r=spec$vectors[2])

spp.scrs <- cbind(spp.scrs, Taxa=tax$Genus)
spp.scrs<-subset(spp.scrs,Taxa!="NA")
spp.sig<-subset(spp.scrs,pVal<=0.01)
spp.sig<-subset(spp.sig,r>=0.3)
spp.sig[,c(1,2)]<-spp.sig[,c(1,2)]*50
spp.plot<-spp.sig[c("ASV_384","ASV_6","ASV_13"),]

scrs4=scrs
spp4=spp.plot

#gw36
otu<-as(otu_table(ps3.gw36.VST),"matrix")
otu<-as.data.frame(otu)
tax=as(tax_table(ps3.gw36.VST),"matrix")
tax<-as.data.frame(tax)
dat<-as(sample_data(ps3.gw36.VST),"data.frame")

nmds.gw36 <- ordinate(ps3.gw36.VST, method="NMDS", distance="euclidean") 
scrs <- as.data.frame(scores(nmds.gw36, display = "sites"))
scrs <- cbind(scrs,Supplementation = dat$Supplementation, Location=dat$Location)
scrs$Location<-as.character(scrs$Location)
scrs$Location<-gsub("0","Linköping",scrs$Location)
scrs$Location<-gsub("1","Stockholm",scrs$Location)
spec <- envfit(nmds.gw36, otu, perm=999)

spp.scrs <- as.data.frame(scores(spec, display="vectors"))
spp.scrs <- cbind(spp.scrs, Species =rownames(spp.scrs))
spp.scrs <- cbind(spp.scrs, pVal=spec$vectors$pvals)
spp.scrs <- cbind(spp.scrs, r=spec$vectors[2])

spp.scrs <- cbind(spp.scrs, Taxa=tax$Genus)
spp.scrs<-subset(spp.scrs,Taxa!="NA")
spp.sig<-subset(spp.scrs,pVal<=0.01)
spp.sig<-subset(spp.sig,r>=0.3)
spp.sig[,c(1,2)]<-spp.sig[,c(1,2)]*50
spp.plot<-spp.sig[c("ASV_384","ASV_6","ASV_13"),]

scrs5=scrs

#2y
otu<-as(otu_table(ps3.2y.VST),"matrix")
otu<-as.data.frame(otu)

tax=as(tax_table(ps3.2y.VST),"matrix")
tax<-as.data.frame(tax)

dat<-as(sample_data(ps3.2y.VST),"data.frame")

nmds.2y <- ordinate(ps3.2y.VST, method="NMDS", distance="euclidean") 
scrs <- as.data.frame(scores(nmds.2y, display = "sites"))
scrs <- cbind(scrs,Supplementation = dat$Supplementation, Location=dat$Location)
scrs$Location<-as.character(scrs$Location)
scrs$Location<-gsub("0","Linköping",scrs$Location)
scrs$Location<-gsub("1","Stockholm",scrs$Location)

scrs6=scrs

## plot ##
p1<-ggplot(scrs1) +
 geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Supplementation, shape=Location)) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("NMDS1")+ylab("NMDS2") + scale_color_manual(values = c("#E69F00","#0072B2")) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp1,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")),size=0.5, colour = "grey",inherit_aes=FALSE) +
  geom_text(data = spp1, aes(x = NMDS1-.02, y = NMDS2-.02, label = Taxa),
            size = 3) +
  ggtitle("1w ***")

p1ellipse<-ggplot(scrs1, aes(x = NMDS1, y = NMDS2, colour = Supplementation, shape=Location)) +
 geom_point() +
 stat_ellipse(aes(group = Location), colour = "black") +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("NMDS1")+ylab("NMDS2") + scale_color_manual(values = c("#E69F00","#0072B2")) +
 coord_fixed() +
 ggtitle("1w ***")

p2<-ggplot(scrs2) +
 geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Supplementation, shape=Location)) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("NMDS1")+ylab("NMDS2") + scale_color_manual(values = c("#E69F00","#0072B2")) +
  coord_fixed() +
  geom_segment(data = spp2,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")),size=0.5, colour = "grey",inherit_aes=FALSE) +
  geom_text(data = spp2, aes(x = NMDS1, y = NMDS2, label = Taxa),
            size = 3) +
  ggtitle("2w ***")

p3<-ggplot(scrs3) +
 geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Supplementation, shape=Location)) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("NMDS1")+ylab("NMDS2") + scale_color_manual(values = c("#E69F00","#0072B2")) +
  coord_fixed() + 
  geom_segment(data = spp3,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")), size=0.5,colour = "grey",inherit_aes=FALSE) +
  geom_text(data = spp3, aes(x = NMDS1-.02, y = NMDS2-.02, label = Taxa),
            size = 3) +
  ggtitle("3w ***")

p4<-ggplot(scrs4) +
 geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Supplementation, shape=Location)) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("NMDS1")+ylab("NMDS2") + scale_color_manual(values = c("#E69F00","#0072B2")) +
  coord_fixed() + 
  geom_segment(data = spp4,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")), size=0.5,colour = "grey",inherit_aes=FALSE) +
  geom_text(data = spp4, aes(x = NMDS1-.02, y = NMDS2-.02, label = Taxa),
            size = 3) +
  ggtitle("4w ***")

p5<-ggplot(scrs5) +
 geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Supplementation, shape=Location)) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("NMDS1")+ylab("NMDS2") + scale_color_manual(values = c("#E69F00","#0072B2")) +
 coord_fixed() +
 ggtitle("PMW36")

p6<-ggplot(scrs6) +
 geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Supplementation, shape=Location)) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("NMDS1")+ylab("NMDS2") + scale_color_manual(values = c("#E69F00","#0072B2")) +
 coord_fixed() +
 ggtitle("2y")

ggarrange(p1+ rremove("legend"), 
  p2+ rremove("legend"), 
  p3+ rremove("legend"), 
  p4+ rremove("legend"), 
  p5+ rremove("legend"), 
  p6+ rremove("legend"),
  align = c("hv"),
  common.legend = TRUE,
  legend="bottom",
  labels="AUTO",
    ncol = 3, nrow = 2)

# stress nmds
stress.validation<-function(x){
  asv=as(otu_table(x), "matrix")
  mds1<-metaMDS(asv,distance="euclidean",k=1,autotransform=F,trymax=200)
  mds2<-metaMDS(asv,distance="euclidean",k=2,autotransform=F,trymax=200)
  mds3<-metaMDS(asv,distance="euclidean",k=3,autotransform=F,trymax=200)
  mds4<-metaMDS(asv,distance="euclidean",k=4,autotransform=F,trymax=200)
  mds5<-metaMDS(asv,distance="euclidean",k=5,autotransform=F,trymax=200)
  mds6<-metaMDS(asv,distance="euclidean",k=6,autotransform=F,trymax=200)
  mds7<-metaMDS(asv,distance="euclidean",k=7,autotransform=F,trymax=200)
  mds8<-metaMDS(asv,distance="euclidean",k=8,autotransform=F,trymax=200)
  mds9<-metaMDS(asv,distance="euclidean",k=9,autotransform=F,trymax=200)
  mds10<-metaMDS(asv,distance="euclidean",k=10,autotransform=F,trymax=200)
  return(stress.val<-c(mds1$stress,mds2$stress,mds3$stress,mds4$stress,mds5$stress,mds6$stress,mds7$stress,mds8$stress,mds9$stress,mds10$stress))
}

# takes long time to run
#s1<-stress.validation(ps3.1w) 
#s2<-stress.validation(ps3.2w)
#s3<-stress.validation(ps3.3w)
#s4<-stress.validation(ps3.4w)
#s5<-stress.validation(ps3.gw36)
#s6<-stress.validation(ps3.2y)

par(mfrow=c(2,3))
plot(1:10,s1,type='b')
plot(1:10,s2,type='b')
plot(1:10,s3,type='b')
plot(1:10,s4,type='b')
plot(1:10,s5,type='b')
plot(1:10,s6,type='b')

# Assess goodness of ordination fit (stress plot)
stressplot(nmds.1w)
stressplot(nmds.2w)
stressplot(nmds.3w)
stressplot(nmds.4w)
stressplot(nmds.gw36)
stressplot(nmds.2y)
stressplot(nmds.1w0)
stressplot(nmds.1w1)

par(mfrow=c(2,2))
plot(1:10,s10,type='b')
plot(1:10,s11,type='b') 
stressplot(nmds.1w0)
stressplot(nmds.1w1)

########################################################################
# ANOSIM
########################################################################
set.seed(60)
anosim_supplementation<-function(x){
  sup = get_variable(x,"Supplementation")
  asv = as(otu_table(x),"matrix")
  ano<-anosim(as.data.frame(asv),grouping=sup,distance="euclidean",permutations=999)
  res<-data.frame(ano$signif,ano$statistic)
  return(res)
}

ano1 = anosim_supplementation(ps3.1w.VST)
ano2 = anosim_supplementation(ps3.2w.VST)
ano3 = anosim_supplementation(ps3.3w.VST)
ano4 = anosim_supplementation(ps3.4w.VST)
ano5 = anosim_supplementation(ps3.gw36.VST)
ano6 = anosim_supplementation(ps3.2y.VST)

ano_sup<-rbind(ano1,ano2,ano3,ano4,ano5,ano6)
rownames(ano_sup)<-c("w1","w2","w3","w4","gw36","y2")
round(ano_sup,2)
#summary(ano1)
#plot(ano1)

# False discovery rate
vec.ano<-as.vector(ano_sup$ano.signif)
p.adjust(vec.ano, method = "BH", n = length(vec.ano))
p.adjust(vec.ano, method = "bonferroni", n = length(vec.ano))

ano_sup_ASV6<-rbind(ano1,ano2,ano3,ano4,ano5)
rownames(ano_sup_ASV6)<-c("w1","w2","w3","w4","gw36")

anosim_location<-function(x){
  sup = get_variable(x,"Location")
  asv = as(otu_table(x),"matrix")
  ano<-anosim(as.data.frame(asv),grouping=sup,distance="euclidean",permutations=999)
  res<-data.frame(ano$signif,ano$statistic)
  return(res)
}

ano1 = anosim_location(ps3.1w.VST)
ano2 = anosim_location(ps3.2w.VST)
ano3 = anosim_location(ps3.3w.VST)
ano4 = anosim_location(ps3.4w.VST)
ano5 = anosim_location(ps3.gw36.VST)
ano6 = anosim_location(ps3.2y.VST)

ano_loc<-rbind(ano1,ano2,ano3,ano4,ano5,ano6)
rownames(ano_loc)<-c("w1","w2","w3","w4","gw36","y2")
#summary(ano1)
#plot(ano1)

#stratified data
ano.1w0 = anosim_supplementation(ps3.1w0.VST) # linköping
ano.1w1 = anosim_supplementation(ps3.1w1.VST) # stockholm

#stratified data
ano.4w0 = anosim_supplementation(ps3.4w0.VST) # AByes
ano.4w1 = anosim_supplementation(ps3.4w1.VST) # ABno

########################################################################
# confounding factors
########################################################################

anosim_gender<-function(x){
  sup = get_variable(x,"Gender")
  asv = as(otu_table(x),"matrix")
  ano<-anosim(as.data.frame(asv),grouping=sup,distance="euclidean",permutations=999)
  res<-data.frame(ano$signif,ano$statistic)
  return(res)
}

g1 = anosim_gender(ps3.1w.VST)
g3 = anosim_gender(ps3.3w.VST)

anosim_delivery<-function(x){
  sup = get_variable(x,"Sectio")
  asv = as(otu_table(x),"matrix")
  ano<-anosim(as.data.frame(asv),grouping=sup,distance="euclidean",permutations=999)
  res<-data.frame(ano$signif,ano$statistic)
  return(res)
}

s3 = anosim_delivery(ps3.3w.VST)

anosim_ABv4<-function(x){
  sup = get_variable(x,"ABv4.0d")
  asv = as(otu_table(x),"matrix")
  ano<-anosim(as.data.frame(asv),grouping=sup,distance="euclidean",permutations=999)
  res<-data.frame(ano$signif,ano$statistic)
  return(res)
}

ab4 = anosim_ABv4(ps3.4w.VST)
