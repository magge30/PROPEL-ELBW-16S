#################################################################################################################################
### FIG 5A: L. reuteri prevalence in L. reuteri and placebo group by time point #################################################
#################################################################################################################################

########## IMPORT & PREPARE DATA ##########

#load data file
data <- read.csv("metadata_github.csv",header=T,sep=",")

# add column with dichotomised qPCR results (0 = negative for L. reuteri, 1 = positive for L. reuteri)
data$qPCR.dich <- data$qPCR
data[data$qPCR.dich!=0 & !is.na(data$qPCR.dich), "qPCR.dich"] <- 1
data$qPCR.dich <- as.factor(data$qPCR.dich) # set the new column as factor


########## L. REUTERI PREVALENCE - RETRIEVE FREQUENCIES FOR INFANTS WITH L. REUTERI-POSITIVE FECES ##########

# function to retrieve frequencies for infants with L. reuteri-positive feces
get.freq <- function(input_data){
  
a <- table(input_data$Supplementation, input_data$qPCR.dich)
ma <- matrix(a, nrow=2)
colnames(ma) <- c("Lreuteri.neg", "Lreuteri.pos")
rownames(ma) <- c("L.reuteri", "Placebo")
ma
freq.1 <- ma[1,1]/sum(ma[1,1],ma[1,2])
freq.2 <- ma[1,2]/sum(ma[1,1],ma[1,2])
freq.3 <- ma[2,1]/sum(ma[2,1],ma[2,2])
freq.4 <- ma[2,2]/sum(ma[2,1],ma[2,2])
freq.lrpos <- c(freq.1, freq.2, freq.3, freq.4)

n.lrpos <- c(ma[1,1], ma[1,2], ma[2,1], ma[2,2])
n.total <- c(sum(ma[1,1],ma[1,2]), sum(ma[1,1],ma[1,2]), sum(ma[2,1],ma[2,2]), sum(ma[2,1],ma[2,2]))

Supplementation <- c("L.reuteri", "L.reuteri", "Placebo", "Placebo")
qPCR.dich <- c(0, 1, 0, 1)

freq.d <- data.frame(Supplementation, qPCR.dich, n.lrpos, n.total, freq.lrpos)
freq.d

return(freq.d)
}

# create tables with frequencies per time point
freq.1w <- get.freq(data[data$Timepoint=="1w",])
freq.1w$Timepoint <- "1w"

freq.2w <- get.freq(data[data$Timepoint=="2w",])
freq.2w$Timepoint <- "2w"

freq.3w <- get.freq(data[data$Timepoint=="3w",])
freq.3w$Timepoint <- "3w"

freq.4w <- get.freq(data[data$Timepoint=="4w",])
freq.4w$Timepoint <- "4w"

freq.PMW36 <- get.freq(data[data$Timepoint=="PMW36",])
freq.PMW36$Timepoint <- "PMW36"

# combine all frequency tables into one data frame
freq.all <- rbind(freq.1w, freq.2w, freq.3w, freq.4w, freq.PMW36)


########## L. REUTERI PREVALENCE - FISHER'S EXACT TEST TO COMPARE L. REUTERI PREVALENCE BETWEEN SUPPLEMENTATION GROUPS ##########

# function for Fisher's exact tests to compare L. reuteri prevalence between supplementation groups per time point
do.fisher.test <- function(input_data){
  
  a <- table(input_data$Supplementation, input_data$qPCR.dich)
  ma <- matrix(a, nrow=2)
  colnames(ma) <- c("Lreut.neg", "Lreut.pos")
  rownames(ma) <- c("L.reuteri", "Placebo")
  ma
  
  fisher.pvals <- c()
  fisher.results <- fisher.test(ma)
  fisher.pvals <- c(fisher.pvals, fisher.results$p.value)
  print(fisher.pvals)
}

# perform Fisher's exact tests per time point
fisher.1w <- do.fisher.test(data[data$Timepoint=="1w",])
fisher.2w <- do.fisher.test(data[data$Timepoint=="2w",])
fisher.3w <- do.fisher.test(data[data$Timepoint=="3w",])
fisher.4w <- do.fisher.test(data[data$Timepoint=="4w",])
fisher.PMW36 <- do.fisher.test(data[data$Timepoint=="PMW36",])

# combine all results (p values) in one vector
fisher.pvalue <- c(fisher.1w, fisher.2w, fisher.3w, fisher.4w, fisher.PMW36)


########## L. REUTERI PREVALENCE - OVERVIEW TABLE ##########

# create nice overview table with L. reuteri prevalence per supplementation group and time point with results from Fisher's exact tests
timepoint <- c("1w","2w","3w","4w","PMW36")

placebo <- freq.all[freq.all$Supplementation=="Placebo" & freq.all$qPCR.dich==1,]
colnames(placebo)[3:5] <- c("placebo.n.lrpos", "placebo.n.total", "placebo.freq.lrpos")

lreuteri <- freq.all[freq.all$Supplementation=="L.reuteri" & freq.all$qPCR.dich==1,]
colnames(lreuteri)[3:5] <- c("lreuteri.n.lrpos", "lreuteri.n.total", "lreuteri.freq.lrpos")

nice.table <- data.frame(timepoint, placebo[,3:5], lreuteri[,3:5], fisher.pvalue)
nice.table$BH.adj.pvalue <- p.adjust(nice.table$fisher.pvalue, method="BH") # adjust p values with the method from Benjamini & Hochberg
nice.table


########## L. REUTERI PREVALENCE - GRAPH FOR FIGURE 5A ##########

### graph showing L. reuteri prevalence by supplementation group and time point

# prepare data for graph
fl.1w <- freq.1w[freq.1w$qPCR.dich==1,]
fl.1w$sup.n <- paste0(fl.1w$Supplementation, " = ", fl.1w$n.total)
fl.1w$sup.n <- gsub("L.reuteri", "Lr", fl.1w$sup.n)
fl.1w$sup.n <- gsub("Placebo", "Pl", fl.1w$sup.n)
fl.1w$sup.n <- ordered(fl.1w$sup.n, levels = c("Pl = 54", "Lr = 54"))

fl.2w <- freq.2w[freq.2w$qPCR.dich==1,]
fl.2w$sup.n <- paste0(fl.2w$Supplementation, "=", fl.2w$n.total)
fl.2w$sup.n <- gsub("L.reuteri", "Lr", fl.2w$sup.n)
fl.2w$sup.n <- gsub("Placebo", "Pl", fl.2w$sup.n)
fl.2w$sup.n <- ordered(fl.2w$sup.n, levels = c("Pl=55", "Lr=54"))

fl.3w <- freq.3w[freq.3w$qPCR.dich==1,]
fl.3w$sup.n <- paste0(fl.3w$Supplementation, "=", fl.3w$n.total)
fl.3w$sup.n <- gsub("L.reuteri", "Lr", fl.3w$sup.n)
fl.3w$sup.n <- gsub("Placebo", "Pl", fl.3w$sup.n)
fl.3w$sup.n <- ordered(fl.3w$sup.n, levels = c("Pl=51", "Lr=51"))

fl.4w <- freq.4w[freq.4w$qPCR.dich==1,]
fl.4w$sup.n <- paste0(fl.4w$Supplementation, "=", fl.4w$n.total)
fl.4w$sup.n <- gsub("L.reuteri", "Lr", fl.4w$sup.n)
fl.4w$sup.n <- gsub("Placebo", "Pl", fl.4w$sup.n)
fl.4w$sup.n <- ordered(fl.4w$sup.n, levels = c("Pl=48", "Lr=53"))

fl.PMW36 <- freq.PMW36[freq.PMW36$qPCR.dich==1,]
fl.PMW36$sup.n <- paste0(fl.PMW36$Supplementation, "=", fl.PMW36$n.total)
fl.PMW36$sup.n <- gsub("L.reuteri", "Lr", fl.PMW36$sup.n)
fl.PMW36$sup.n <- gsub("Placebo", "Pl", fl.PMW36$sup.n)
fl.PMW36$sup.n <- ordered(fl.PMW36$sup.n, levels = c("Pl=41", "Lr=50"))

fl.all <- rbind(fl.1w, fl.2w, fl.3w, fl.4w, fl.PMW36)
fl.all

# prepare labels for graph, set colors etc.
Supplementation.colors <- c("Placebo" = "#E69F00", "L.reuteri" = "#0072B2")
fl.all$Timepoint <- ordered(fl.all$Timepoint, levels = c("1w", "2w", "3w", "4w", "PMW36"))
label_y<-expression(paste("Percentage of infants with ",italic("L. reuteri"),"-positive stool sample"))
stripnames <- c(paste0("\n 1w \n"),paste0("\n 2w \n"),paste0("\n 3w \n"),paste0("\n 4w \n"),paste0("\n PMW36 \n"))
names(stripnames) <- c("1w", "2w", "3w", "4w", "PMW36")
annotation <- data.frame(sup.n = c(rep(1.5,5)),
                         freq.lrpos = c(rep(1.05,5)),
                         Supplementation = c("Placebo","Placebo","Placebo","Placebo","Placebo"),
                         label_x = c("***","***","***","***","***"),
                         Timepoint = c("1w", "2w", "3w", "4w", "PMW36"))
label_x <- c("***","***","***","***","***")

# graph
library(ggplot2)
g1 <- ggplot(fl.all,
       aes(x=sup.n, y=freq.lrpos, fill=Supplementation))+
  geom_bar(stat = "identity",color="black",
           position = position_dodge(width = 1))+
  facet_grid(.~Timepoint, scales = "free_x", labeller=labeller(Timepoint=stripnames))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="", y=label_y)+
  theme(legend.position = "none")+
  scale_fill_manual(values=Supplementation.colors)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),
                     labels = c("0%","20%","40%","60%","80%","100%"))+
  theme(strip.text.x=element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20))+
  geom_text(data=annotation, label=label_x, size=7)
g1


#################################################################################################################################
### FIG 5B: L. reuteri abundance in L. reuteri and placebo group by time point #################################################
#################################################################################################################################

########## IMPORT & PREPARE DATA ##########

#load data file
data <- read.csv("metadata_github.csv",header=T,sep=",")


########## CHECK DATA DISTRIBUTIONS ########## 

library(ggplot2)

# function to create histograms and qqplots for L. reuteri abundance data per group and time point
distribution.check <- function(input_data, output_name){
  pdf(paste0("histograms_", output_name, ".pdf"))
  for(i in 6){
    #make histograms to see data distributions and save them in pdf
    histograms<-ggplot(input_data,aes(x=as.numeric(input_data[,i])))+
      geom_histogram()+
      theme_classic()+
      ggtitle("L.reuteri.abundance")
    print(histograms)
  }
  dev.off()
  
  pdf(paste0("qqplots_", output_name, ".pdf"))
  for(i in 6){
    qqplots<-ggplot(input_data,aes(sample=as.numeric(input_data[,i])))+
      stat_qq()+
      stat_qq_line()+
      theme_classic()+
      ggtitle("L.reuteri.abundance")
    print(qqplots)
  }
  dev.off()
}

## checking levels in lr-negative and lr-positive infants
distribution.check(data[data$Timepoint== "1w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),], "1w.placebo")
distribution.check(data[data$Timepoint== "1w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),], "1w.lreuteri")

distribution.check(data[data$Timepoint== "2w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),], "2w.placebo")
distribution.check(data[data$Timepoint== "2w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),], "2w.lreuteri")

distribution.check(data[data$Timepoint== "3w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),], "3w.placebo")
distribution.check(data[data$Timepoint== "3w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),], "3w.lreuteri")

distribution.check(data[data$Timepoint== "4w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),], "4w.placebo")
distribution.check(data[data$Timepoint== "4w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),], "4w.lreuteri")

distribution.check(data[data$Timepoint== "PMW36" & data$Supplementation=="Placebo" & !is.na(data$qPCR),], "PMW36.placebo")
distribution.check(data[data$Timepoint== "PMW36" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),], "PMW36.lreuteri")

# --> L. reuteri abundance data is left skewed
# --> present data as median (IQR) and use non-parametric statistical tests


########## L. REUTERI ABUNDANCE - RETRIEVE MEDIANS AND INTERQUARTILE RANGES (IQR) ########## 

# function to retreive median and IQR for L. reuteri abundance per supplementation group and time point
get.median.iqr.n <- function(input_data, rowname){
  placebo <- input_data[input_data$Supplementation=="Placebo",]
  lreuteri <- input_data[input_data$Supplementation=="L.reuteri",]
  m.p <- median(placebo[,6])
  m.l <- median(lreuteri[,6])
  iqr.025.p <- quantile(placebo[,6], c(0.25))
  iqr.075.p <- quantile(placebo[,6], c(0.75))
  iqr.025.l <- quantile(lreuteri[,6], c(0.25))
  iqr.075.l <- quantile(lreuteri[,6], c(0.75))
  n.p <- nrow(placebo)
  n.l <- nrow(lreuteri)
  merged <- cbind(m.p, iqr.025.p, iqr.075.p, n.p,
                  m.l, iqr.025.l, iqr.075.l, n.l)
  colnames(merged) <- c("median_p", "iqr_025_p", "iqr_075_p", "n_p",
                        "median_l", "iqr_025_l", "iqr_075_l", "n_l")
  rownames(merged) <- paste0(rowname)
  return(merged)
}

# retrieve medians and IQR per time point
med.iqr.1w <- get.median.iqr.n(data[data$Timepoint=="1w" & !is.na(data$qPCR),], "1w")
med.iqr.2w <- get.median.iqr.n(data[data$Timepoint=="2w" & !is.na(data$qPCR),], "2w")
med.iqr.3w <- get.median.iqr.n(data[data$Timepoint=="3w" & !is.na(data$qPCR),], "3w")
med.iqr.4w <- get.median.iqr.n(data[data$Timepoint=="4w" & !is.na(data$qPCR),], "4w")
med.iqr.PMW36 <- get.median.iqr.n(data[data$Timepoint=="PMW36" & !is.na(data$qPCR),], "PMW36")

# combine all medians and IQR into one data frame
med.iqr.all <- data.frame(rbind(med.iqr.1w, med.iqr.2w, med.iqr.3w, med.iqr.4w, med.iqr.PMW36))


########## L. REUTERI ABUNDANCE - MANN-WHITNEY U TESTS TO COMPARE L. REUTERI ABUNDANCE BETWEEN SUPPLEMENTATION GROUP ########## 

# function for Mann-Whitney U test
mwu.test <- function(pl,lr){
  mwu.pvals <- c()
  mwu.results <- wilcox.test(pl, lr)
  mwu.pvals <- c(mwu.pvals, mwu.results$p.value)
  print(mwu.pvals)
}

# perform Mann-Whitney U tests per time point
mwu.1w <- mwu.test(data[data$Timepoint=="1w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),6], data[data$Timepoint=="1w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),6])
mwu.2w <- mwu.test(data[data$Timepoint=="2w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),6], data[data$Timepoint=="2w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),6])
mwu.3w <- mwu.test(data[data$Timepoint=="3w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),6], data[data$Timepoint=="3w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),6])
mwu.4w <- mwu.test(data[data$Timepoint=="4w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),6], data[data$Timepoint=="4w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),6])
mwu.PMW36 <- mwu.test(data[data$Timepoint=="PMW36" & data$Supplementation=="Placebo" & !is.na(data$qPCR),6], data[data$Timepoint=="PMW36" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),6])

# combine all results (p values) in one vector
mwu_pvalue <- c(mwu.1w, mwu.2w, mwu.3w, mwu.4w, mwu.PMW36)


########## L. REUTERI ABUNDANCE - OVERVIEW TABLE ##########

# create nice overview table with L. reuteri abundance per supplementation group and time point with results from Mann-Whitney U tests
med.iqr.all$mwu_pvalue <- mwu_pvalue # add colum with p values from Mann-Whitney U tests to data frame with medians and IQR
med.iqr.all$BH_adj_pvalue <- p.adjust(med.iqr.all$mwu_pvalue, method="BH") # adjust p values with the method from Benjamini & Hochberg
med.iqr.all


########## L. REUTERI ABUNDANCE - KRUSKAL-WALLIS TESTS TO COMPARE L. REUTERI ABUNDANCE IN L. REUTERI GROUP OVER TIME (NEONATAL TIME POINTS ONLY) ########## 

# perform Kruskal-Wallis test
library(PMCMRplus)
kw.results <- kruskalTest(qPCR ~ Timepoint, data=data[data$Supplementation=="L.reuteri" & data$Timepoint!="2y",])
kw.results$p.value
# --> p-value = 0.000008222
# --> significant, perform posthoc test

# Dunn's posthoc test with Benjamini-Hochberg-adjusted p value
library(FSA)
dt.results <- dunnTest(qPCR ~ Timepoint, data=data[data$Supplementation=="L.reuteri" & data$Timepoint!="2y",], method="bh")
dt.results <- dt.results$res
dt.results[c(1,2,4,7,3,5,8,6,9,10),] # show results from Dunn posthoc test sorted by time
# --> significant differences for comparing L. reuteri abundances between 1w-PMW36, 2w-PMW36, 3w-PMW36 and 4w-PMW36


########## L. REUTERI ABUNDANCE - GRAPH FOR FIGURE 5B ########## 

### graph showing L. reuteri abundance by supplementation group and time point

# prepare data for graph
pl.1w <- data[data$Timepoint=="1w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),]
lr.1w <- data[data$Timepoint=="1w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),]
pl.1w$sup.n <- as.factor(paste0("Pl=", nrow(pl.1w)))
lr.1w$sup.n <- as.factor(paste0("Lr=", nrow(lr.1w)))
pl.lr.1w <- rbind(pl.1w, lr.1w)
pl.lr.1w$sup.n <- ordered(pl.lr.1w$sup.n, levels = c("Pl=54", "Lr=54"))

pl.2w <- data[data$Timepoint=="2w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),]
lr.2w <- data[data$Timepoint=="2w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),]
pl.2w$sup.n <- as.factor(paste0("Pl=", nrow(pl.2w)))
lr.2w$sup.n <- as.factor(paste0("Lr-", nrow(lr.2w)))
pl.lr.2w <- rbind(pl.2w, lr.2w)
pl.lr.2w$sup.n <- ordered(pl.lr.2w$sup.n, levels = c("Pl=55", "Lr-54"))

pl.3w <- data[data$Timepoint=="3w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),]
lr.3w <- data[data$Timepoint=="3w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),]
pl.3w$sup.n <- as.factor(paste0("Pl=", nrow(pl.3w)))
lr.3w$sup.n <- as.factor(paste0("Lr=", nrow(lr.3w)))
pl.lr.3w <- rbind(pl.3w, lr.3w)
pl.lr.3w$sup.n <- ordered(pl.lr.3w$sup.n, levels = c("Pl=51", "Lr=51"))

pl.4w <- data[data$Timepoint=="4w" & data$Supplementation=="Placebo" & !is.na(data$qPCR),]
lr.4w <- data[data$Timepoint=="4w" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),]
pl.4w$sup.n <- as.factor(paste0("Pl=", nrow(pl.4w)))
lr.4w$sup.n <- as.factor(paste0("Lr=", nrow(lr.4w)))
pl.lr.4w <- rbind(pl.4w, lr.4w)
pl.lr.4w$sup.n <- ordered(pl.lr.4w$sup.n, levels = c("Pl=48", "Lr=53"))

pl.PMW36 <- data[data$Timepoint=="PMW36" & data$Supplementation=="Placebo" & !is.na(data$qPCR),]
lr.PMW36 <- data[data$Timepoint=="PMW36" & data$Supplementation=="L.reuteri" & !is.na(data$qPCR),]
pl.PMW36$sup.n <- as.factor(paste0("Pl=", nrow(pl.PMW36)))
lr.PMW36$sup.n <- as.factor(paste0("Lr=", nrow(lr.PMW36)))
pl.lr.PMW36 <- rbind(pl.PMW36, lr.PMW36)
pl.lr.PMW36$sup.n <- ordered(pl.lr.PMW36$sup.n, levels = c("Pl=41", "Lr=50"))

pl.lr.all <- rbind(pl.lr.1w, pl.lr.2w, pl.lr.3w, pl.lr.4w, pl.lr.PMW36)


# set all infants with a fecal sample negative for L. reuteri to '1' instead of '0' so that they can be displayed on a log scale
plot.data <- pl.lr.all
plot.data[plot.data$qPCR==0,6] <- 1

# prepare labels for graph, set colors etc.
options(scipen = 999)
Supplementation.colors <- c("Placebo" = "#D55E00", "L.reuteri" = "#56B4E9")
Supplementation.fill <- c("Placebo" = "#E69F00", "L.reuteri" = "#0072B2")
plot.data$Supplementation <- ordered(plot.data$Supplementation, c("Placebo", "L.reuteri"))
plot.data$Timepoint <- ordered(plot.data$Timepoint, c("1w", "2w", "3w", "4w", "PMW36"))
label_y <- expression(paste(italic("L. reuteri"), " bacteria per 1 g wet feces"))
stripnames <- c(paste0("\n 1w \n"),paste0("\n 2w \n"),paste0("\n 3w \n"),paste0("\n 4w \n"),paste0("\n PMW36 \n"))
names(stripnames) <- c("1w", "2w", "3w", "4w", "PMW36")

# graph
library(ggplot2)
g1 <- ggplot(plot.data, aes(x=sup.n, y=qPCR, fill=Supplementation))+
  geom_boxplot(position = position_dodge(1, preserve="single"), outlier.size=2)+
  geom_point(position="jitter", aes(color=Supplementation, shape="."))+
  facet_grid(.~Timepoint, scales="free_x", labeller=labeller(Timepoint=stripnames))+
  theme_bw()+
  labs(x="", y=label_y)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size=20),
        panel.grid.minor = element_blank())+
  scale_y_log10(limits = c(1e-1, 8e+9),
                breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8, 1e+9),
                labels = c("10^0", "10^1", "10^2", "10^3", "10^4", "10^5", "10^6", "10^7", "10^8", "10^9"))+
  scale_fill_manual(values=Supplementation.fill)+
  scale_color_manual(values=Supplementation.colors)
g1





















