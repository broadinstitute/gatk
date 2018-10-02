#!/usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

# java -jar /seq/software/picard/current/bin/picard-private.jar UpdateVcfSequenceDictionary INPUT=/dsde/data/deep/vqsr/vcfs/nist_na12878_giab_hg38.vcf.gz OUTPUT=/dsde/data/deep/vqsr/vcfs/nist_na12878_giab_hg38_sd_fix.vcf SD=/dsde/data/deep/vqsr/vcfs/illumina_na12878_platinum_scored_chr2.vcf
# ./gatk GenotypeConcordance --CALL_VCF=/dsde/data/deep/vqsr/vcfs/illumina_na12878_platinum_scored_chr20.vcf.gz --CALL_SAMPLE=NA12878 --TRUTH_VCF=/dsde/data/deep/vqsr/vcfs/nist_na12878_giab_hg38_sd_fix.vcf.gz --TRUTH_SAMPLE=HG001 --INTERVALS=/dsde/data/deep/vqsr/beds/HG001_NA12878_GRCh38_GIAB_highconf.interval_list --OUTPUT_VCF=true -O=na12878_hg38_chr20_concordance
# ./gatk VariantsToTable -V /dsde/data/deep/vqsr/vcfs/illumina_na12878_platinum_scored_chr2.vcf.gz -F CHROM -F POS -F REF -F ALT -F FILTER -F G947_SITE_LABELLED_RRAB -F EVENTLENGTH -F AC -F MULTI-ALLELIC -F TRANSITION -F TYPE --show-filtered -O ~/Documents/illumin_chr2.table
# ./gatk VariantsToTable -V na12878_hg38_chr2_concordance.genotype_concordance.vcf.gz -F CHROM -F POS -F REF -F ALT -F CONC_ST -O ~/Documents/illumina_chr2_truth.table
# d <- read.table("illumina_chr2.table", header=TRUE)
# dt <- read.table("illumina_chr2_truth.table", header=TRUE)
#score_key <- "G947_SITE_LABELLED_RRAB"
#d <- read.table("g947h_chr20.table", header=TRUE)
#dt <- read.table("g947_chr20_truth.table", header=TRUE)
# d <- read.table("g94982_chr20.table", header=TRUE)
# dt <- read.table("g94982_chr20_truth.table", header=TRUE)
# score_key <- "CNN_2D"
# plot_title <- "yo"

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("We need 3 arguments: call_vcf_table concordance_vcf_table score_key")
}

print("try to load VCF table.")
d <- read.table(args[1], header=TRUE)
print("try to load VCF Truth table.")
dt <- read.table(args[2], header=TRUE)
score_key <- args[3]
score_label <- paste(score_key, " LOD Score")
plot_title <- gsub( ".vcf.gz.table", "", basename(args[1]))
num_bins <- 50
bin_by_quantile <-  FALSE

get_proportion <- function(d, num_bins, column_to_sum, quality_column) {
  x <- rowsum(column_to_sum, quality_column, na.rm =T)
  idx <- row.names(x)
  
  for (i in 1:num_bins) {
    qsum <- sum(quality_column==as.numeric(idx[i]))
    if (!is.na(x[i]) && qsum>0) {
      x[i] <- x[i] / qsum
    }
  }
  return(x[quality_column])
}

print("try to merge.")
d <- merge(d, dt, by=c("CHROM", "POS", "REF", "ALT"))
d$TP <- as.numeric(d$CONC_ST!="FP,TN" & d$CONC_ST!="FP" & d$CONC_ST!="EMPTY")
d$True_Positive <- d$CONC_ST!="FP,TN" & d$CONC_ST!="FP" & d$CONC_ST!="EMPTY"
d$Unfiltered <- d$FILTER == "PASS" | d$FILTER == "."
d$SNP <- d$EVENTLENGTH == 0
d$ONE <- 1
x <- rowsum(d$ONE, d$EVENTLENGTH)
d$EVENTLENGTH_SUM <- x[as.factor(d$EVENTLENGTH)]
d$Variant_Type <- paste(d$TYPE, as.factor(d$EVENTLENGTH<0))
d$Truth_Status <- ifelse(d$True_Positive & d$Unfiltered, "True Positive", ifelse(d$True_Positive & !d$Unfiltered, "False Negative", ifelse(!d$True_Positive & d$Unfiltered, "False Positive", "True Negative")))
statusColor <- c("True Positive" = "springgreen3", "True Negative" = "aquamarine4", "False Positive" = "red", "False Negative" = "orange")

# All variant plots
print("Make all variant plots.")
p1 <- ggplot(d, aes(get(score_key), color=SNP, fill=SNP)) + 
  scale_fill_discrete(name="Variant\nType", breaks=c("TRUE", "FALSE"), labels=c("SNPs", "INDELs")) +
  geom_density(alpha=0.55) +
  ggtitle(plot_title) +
  xlab(score_label) +
  guides(color=FALSE)

p2 <- ggplot(d, aes(x=EVENTLENGTH, y=get(score_key), group=EVENTLENGTH, color=Truth_Status, shape=Variant_Type)) + 
  scale_color_manual(values=statusColor) +
  scale_shape_discrete(name='', breaks=c("INDEL TRUE", "INDEL FALSE", "SNP FALSE"), labels=c("Deletion", "Insertion", "SNP")) +
  geom_jitter(height = 0, width = 0.1, alpha=0.6) + 
  ggtitle(plot_title) +
  ylab(score_label) + 
  xlab("Event Length: - Deletions, 0 SNPs, + Insertions")

p3 <- ggplot(d, aes(x=EVENTLENGTH, y=get(score_key), group=EVENTLENGTH, color=Truth_Status)) + xlim(-20, 20) + 
  scale_color_manual(values=statusColor) +
  geom_jitter(height = 0, width = 0.1, alpha=0.4) + 
  geom_violin(color="grey", alpha=0) + 
  geom_text(aes(x=EVENTLENGTH, y=14, label=EVENTLENGTH_SUM), color="grey30", size=2, angle=60) +
  ggtitle(plot_title) +
  ylab(score_label) + 
  xlab("Event Length: - Deletions, 0 SNPs, + Insertions")

p4 <- ggplot(d, aes(x=EVENTLENGTH, y=get(score_key), group=EVENTLENGTH, color=Truth_Status)) + xlim(-10, 10) + 
  scale_color_manual(values=statusColor) +
  geom_jitter(height = 0, width = 0.2, alpha=0.4) + 
  geom_violin(color="grey", alpha=0) + 
  geom_text(aes(x=EVENTLENGTH, y=14, label=EVENTLENGTH_SUM), color="grey30", size=3, angle=30) +
  ylab(score_label) + 
  xlab("Event Length: - Deletions, 0 SNPs, + Insertions")

p5 <- ggplot(d, aes(x=EVENTLENGTH, y=get(score_key), group=EVENTLENGTH, color=Truth_Status)) + 
  scale_color_manual(values=statusColor) +
  xlim(-5, 5) + 
  geom_jitter(height = 0, width = 0.35, alpha=0.4) + 
  geom_violin(color="grey", alpha=0.0) + 
  geom_text(aes(x=EVENTLENGTH, y=14, label=EVENTLENGTH_SUM), color="grey30", size=4, angle=30) +
  ggtitle(plot_title) +
  ylab(score_label) + 
  xlab("Event Length: - Deletions, 0 SNPs, + Insertions")

# SNP specific plots
print("Make SNP plots.")
snps <- subset(d, EVENTLENGTH == 0)
my_breaks <- ifelse(bin_by_quantile, quantile(snps[[score_key]], probs = seq(0, 1, 1.0/num_bins), na.rm=T), num_bins)
snps$QUALITY_BIN <- cut(snps[[score_key]], breaks=my_breaks, include.lowest=T, labels=F)
snps$QUALITY_BIN_RANGE <- cut(snps[[score_key]], breaks=my_breaks, include.lowest=T)
mine <- lapply(strsplit(sapply(levels(snps$QUALITY_BIN_RANGE), function(x) substr(x, 2, nchar(x)-1)), ","), as.numeric)
df <- data.frame(matrix(unlist(mine), nrow=num_bins, byrow=T))
q_means <- rowMeans(df)
snps$QUALITY_LOD <- q_means[snps$QUALITY_BIN] 
snps$TPR_PREDICTION <- exp(snps$QUALITY_LOD) / (1 + exp(snps$QUALITY_LOD) )

x <- rowsum(snps$ONE, snps$QUALITY_BIN)
snps$BIN_SUM <- x[snps$QUALITY_BIN]
snps$TRANSVERSION <- as.numeric( abs(snps$TRANSITION)==0 )
snps$TPR <- get_proportion(snps, num_bins, snps$TP, snps$QUALITY_BIN)
ti <- get_proportion(snps, num_bins, snps$TRANSITION, snps$QUALITY_BIN)
tv <- get_proportion(snps, num_bins, snps$TRANSVERSION, snps$QUALITY_BIN)
snps$TI_TV <- ti/tv

p6 <- ggplot(snps, aes(x=get(score_key), y=TI_TV, group=QUALITY_BIN, color=Truth_Status, shape=TRANSITION==1)) + 
  scale_color_manual(values=statusColor) +
  scale_shape_discrete(name='', breaks=c("TRUE", "FALSE"), labels=c("Transition", "Transversion")) +
  geom_point() + 
  geom_line(color="grey") + 
  ggtitle("Transition Transversion Ratio per score bin") +
  xlab(score_label) +
  ylim(0, 4)

p7 <- ggplot(snps, aes(x=TPR_PREDICTION, y=TPR, group=QUALITY_BIN, color=Truth_Status)) + 
  scale_color_manual(values=statusColor) +
  geom_jitter(height = 0.01, width = 0.01, alpha=0.4) + 
  ggtitle(paste("SNP Calibration", plot_title)) +
  ylim(0, 1) + xlim(0, 1)

# INDEL specific plots
print("Make INDEL plots.")
indels <- subset(d, EVENTLENGTH != 0)
my_breaks <- ifelse(bin_by_quantile, quantile(indels[[score_key]], probs = seq(0, 1, 1.0/num_bins), na.rm=T), num_bins)
indels$QUALITY_BIN <- cut(indels[[score_key]], breaks=my_breaks, include.lowest=T, labels=F)
indels$QUALITY_BIN_RANGE <- cut(indels[[score_key]], breaks=my_breaks, include.lowest=T)
mine <- lapply(strsplit(sapply(levels(indels$QUALITY_BIN_RANGE), function(x) substr(x, 2, nchar(x)-1)), ","), as.numeric)
df <- data.frame(matrix(unlist(mine), nrow=num_bins, byrow=T))
q_means <- rowMeans(df)
indels$QUALITY_LOD <- q_means[indels$QUALITY_BIN] 
indels$TPR_PREDICTION <- exp(indels$QUALITY_LOD) / (1 + exp(indels$QUALITY_LOD))
x <- rowsum(indels$ONE, indels$QUALITY_BIN)
indels$BIN_SUM <- x[indels$QUALITY_BIN]
indels$TPR <- get_proportion(indels, num_bins, indels$TP, indels$QUALITY_BIN)
indels$ONEBP <- as.numeric(abs(indels$EVENTLENGTH)==1)
indels$PROPORTION_ONEBP <- get_proportion(indels, num_bins, indels$ONEBP, indels$QUALITY_BIN)

p8 <- ggplot(indels, aes(x=get(score_key), y=PROPORTION_ONEBP, group=QUALITY_BIN, color=Truth_Status, shape=EVENTLENGTH<0)) +
  scale_color_manual(values=statusColor) +
  scale_shape_discrete(name='', breaks=c("TRUE", "FALSE"), labels=c("Deletion", "Insertion")) +
  geom_jitter(height = 0.005, width = 0.0, alpha=0.6) +
  geom_line(color="grey") + 
  ggtitle("Proportion of 1bp INDELs per score bin") +
  xlab(score_label)

p9 <- ggplot(indels, aes(x=TPR_PREDICTION, y=TPR, group=QUALITY_BIN, color=Truth_Status)) +
  scale_color_manual(values=statusColor) +
  geom_jitter(height = 0.01, width = 0.01, alpha=0.4) + 
  ggtitle(paste("INDEL Calibration", plot_title)) +
  ylim(0, 1) + xlim(0, 1)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
ggsave(plot=multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9, cols=2), filename = paste(plot_title, "_plots.png", sep=""), width=16, height=22)
