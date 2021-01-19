library("diverse")

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

# extract phyloseq otu and data
otu1<-as(otu_table(ps3.1w),"matrix")
otu2<-as(otu_table(ps3.2w),"matrix")
otu3<-as(otu_table(ps3.3w),"matrix")
otu4<-as(otu_table(ps3.4w),"matrix")
otu5<-as(otu_table(ps3.gw36),"matrix")
otu6<-as(otu_table(ps3.2y),"matrix")

data1<-as(sample_data(ps3.1w),"data.frame")
data2<-as(sample_data(ps3.2w),"data.frame")
data3<-as(sample_data(ps3.3w),"data.frame")
data4<-as(sample_data(ps3.4w),"data.frame")
data36<-as(sample_data(ps3.gw36),"data.frame")
data2y<-as(sample_data(ps3.2y),"data.frame")

# compute alpha diversity 
adt<-function(x){ alfa<-diversity(t(x),type=c("entropy","variety","evenness","berger-parker")); return(alfa) }
alfa1<-adt(otu1)
alfa2<-adt(otu2)
alfa3<-adt(otu3)
alfa4<-adt(otu4)
alfa5<-adt(otu5)
alfa6<-adt(otu6)

# merge diversity index with metada
div1<-merge(alfa1,data1,by="row.names",all.x=TRUE)
div2<-merge(alfa2,data2,by="row.names",all.x=TRUE)
div3<-merge(alfa3,data3,by="row.names",all.x=TRUE)
div4<-merge(alfa4,data4,by="row.names",all.x=TRUE)
div5<-merge(alfa5,data36,by="row.names",all.x=TRUE)
div6<-merge(alfa6,data2y,by="row.names",all.x=TRUE)

# fix data frame
fix<-function(x){
  row.names(x)<-x$"Row.names"
  x = x[,c("entropy","variety","evenness","berger.parker.D","Supplementation","Timepoint")] 
  colnames(x)<-c("Diversity","Richness","Evenness","Dominance","Supplementation","Timepoint")
  return(x)
}

div1<-fix(div1)
div2<-fix(div2)
div3<-fix(div3)
div4<-fix(div4)
div5<-fix(div5)
div6<-fix(div6)

div1m<-melt(div1[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div2m<-melt(div2[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div3m<-melt(div3[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div4m<-melt(div4[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div5m<-melt(div5[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div6m<-melt(div6[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])

# normality test
adt<-function(x){ test<-ad.test(x); return(test$p.value) }
#store the p-values for each column in a separate variable
p.adt1<-apply(div1[,c(1:4)],MARGIN = 2, adt)
p.adt2<-apply(div2[,c(1:4)],MARGIN = 2, adt)
p.adt3<-apply(div3[,c(1:4)],MARGIN = 2, adt)
p.adt4<-apply(div4[,c(1:4)],MARGIN = 2, adt)
p.adt5<-apply(div5[,c(1:4)],MARGIN = 2, adt)
p.adt6<-apply(div6[,c(1:4)],MARGIN = 2, adt)

# mann-whitney
mnt1<-sapply(div1[,c("Diversity","Richness","Evenness","Dominance")], function(x) wilcox.test(x ~ div1$Supplementation)$p.value)
mnt2<-sapply(div2[,c("Diversity","Richness","Evenness","Dominance")], function(x) wilcox.test(x ~ div2$Supplementation)$p.value)
mnt3<-sapply(div3[,c("Diversity","Richness","Evenness","Dominance")], function(x) wilcox.test(x ~ div3$Supplementation)$p.value)
mnt4<-sapply(div4[,c("Diversity","Richness","Evenness","Dominance")], function(x) wilcox.test(x ~ div4$Supplementation)$p.value)
mnt5<-sapply(div5[,c("Diversity","Richness","Evenness","Dominance")], function(x) wilcox.test(x ~ div5$Supplementation)$p.value)
mnt6<-sapply(div6[,c("Diversity","Richness","Evenness","Dominance")], function(x) wilcox.test(x ~ div6$Supplementation)$p.value)

mnt.df<-rbind(mnt1,mnt2,mnt3,mnt4,mnt5,mnt6)
mnt.df<-as.data.frame(mnt.df)

# p.adjust
vec.diversity<-as.vector(mnt.df$Diversity)
vec.richness<-as.vector(mnt.df$Richness)
vec.evenness<-as.vector(mnt.df$Evenness)
all.index<-c(vec.diversity,vec.richness,vec.evenness)
p.adjust(all.index, method = "BH", n = length(all.index))

## Plot ##
library(tidyverse)
library(hrbrthemes)
library(viridis)

cbbPalette<-c("#E69F00","#0072B2")

div1[["Timepoint"]]<-"1w"
div2[["Timepoint"]]<-"2w"
div3[["Timepoint"]]<-"3w"
div4[["Timepoint"]]<-"4w"
div5[["Timepoint"]]<-"PMW36"
div6[["Timepoint"]]<-"2y"

div1m<-melt(div1[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div2m<-melt(div2[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div3m<-melt(div3[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div4m<-melt(div4[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div5m<-melt(div5[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])
div6m<-melt(div6[,c("Diversity","Richness","Evenness","Supplementation","Timepoint")])

# to get first placebo and then L. reuteri in the x-axis
div1m$test <- div1m$Supplementation 
div1m$test <- gsub('Placebo', 'Pl=54', div1m$test)
div1m$test <- gsub('L.reuteri', 'Lr=54', div1m$test)
div1m$test = factor(div1m$test, levels=c("Pl=54","Lr=54")) 

div2m$test <- div2m$Supplementation 
div2m$test <- gsub('Placebo', 'Pl=55', div2m$test)
div2m$test <- gsub('L.reuteri', 'Lr=54.', div2m$test)
div2m$test = factor(div2m$test, levels=c("Pl=55","Lr=54."))

div3m$test <- div3m$Supplementation
div3m$test <- gsub('Placebo', 'Pl=51', div3m$test)
div3m$test <- gsub('L.reuteri', 'Lr=51', div3m$test)
div3m$test = factor(div3m$test, levels=c("Pl=51","Lr=51"))

div4m$test <- div4m$Supplementation
div4m$test <- gsub('Placebo', 'Pl=48', div4m$test)
div4m$test <- gsub('L.reuteri', 'Lr=53', div4m$test)
div4m$test = factor(div4m$test, levels=c("Pl=48","Lr=53"))

div5m$test <- div5m$Supplementation
div5m$test <- gsub('Placebo', 'Pl=41', div5m$test)
div5m$test <- gsub('L.reuteri', 'Lr=50', div5m$test)
div5m$test = factor(div5m$test, levels=c("Pl=41","Lr=50"))

div6m$test <- div6m$Supplementation
div6m$test <- gsub('Placebo', 'Pl=27', div6m$test)
div6m$test <- gsub('L.reuteri', 'Lr=20', div6m$test)
div6m$test = factor(div6m$test, levels=c("Pl=27","Lr=20"))

# merge data
divm<-rbind(div1m,div2m,div3m,div4m,div5m,div6m)

# to sort out facet grid by time point
divm$Timepoint = factor(divm$Timepoint, levels=c("1w","2w","3w","4w","PMW36","2y"))

palpha<-ggplot(divm, aes(x=test, y=value,fill=Supplementation)) + 
  geom_boxplot() + ylab("Alpha Diversity Measure") + theme_bw( ) + facet_grid(variable ~ Timepoint, scales="free",)  +
  scale_fill_manual(values=cbbPalette)  +
  theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.title.x=element_blank(),legend.position="none") 

ggsave("Fig2_alpha.pdf",palpha)




