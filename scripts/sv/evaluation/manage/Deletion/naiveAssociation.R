## This script tries to look at, naively, what type of variables are associated 
##  with the overlapping calls between PacBio and GATK, and what's associated with 
##  those that are present in PacBio's call set but not in GATK's call set, 
##  and vice versa.
library("plotrix")
library("plyr")
library("stringr")

args<-commandArgs(TRUE)
wd <- args[1]
setwd(wd)

################################################################################
# load two background data, manta call set and gatk call set
#   and load intersections
background_pacbio <- read.table("PacBio_primaryContigs_cleanDel.vcf", 
                                stringsAsFactors = F, header = F)
names(background_pacbio) <- c("CHROM", "POS", "ID", "REF", "ALT", 
                              "QUAL", "FILTER", "INFO")

background_manta <- read.table("Manta_PASS_PRECISE_nonBND_primaryContigs_cleanDel.vcf", 
                               stringsAsFactors = F, header = F)
names(background_manta) <- c("CHROM", "POS", "ID", "REF", "ALT", 
                             "QUAL", "FILTER", "INFO", 
                             "FORMAT", "SAMPLE")

background_gatk <- read.table("GATK_primaryContigs_cleanDel.vcf", 
                              stringsAsFactors = F, header = F)
names(background_gatk) <- c("CHROM", "POS", "ID", "REF", "ALT", 
                            "QUAL", "FILTER", "INFO")
################################################################################
gatkIntersectManta <- read.table("Results/cleanDel.gatkVSmanta.intersection_f08r.txt", 
                                 stringsAsFactors = F, header = F, 
                                 na.strings = c(".", "-1", "0"))
names(gatkIntersectManta) <- c("gatkCHR", "gatkPOS", "gatkEND", "gatkAnnot", 
                               "mantaCHR", "mantaPOS", "mantaEND", "mantaAnnot", 
                               "ovpLen")

gatkIntersectPacbio <- read.table("Results/cleanDel.gatkVSpacbio.intersection_f08r.txt", 
                                  stringsAsFactors = F, header = F, 
                                  na.strings = c(".", "-1", "0"))
names(gatkIntersectPacbio) <- c("gatkCHR", "gatkPOS", "gatkEND", "gatkAnnot", 
                                "pbCHR", "pbPOS", "pbEND", "pbAnnot", 
                                "ovpLen")

mantaIntersectPacbio <- read.table("Results/cleanDel.mantaVSpacbio.intersection_f08r.txt", 
                                   stringsAsFactors = F, header = F, 
                                   na.strings = c(".", "-1", "0"))
names(mantaIntersectPacbio) <- c("mantaCHR", "mantaPOS", "mantaEND", "mantaAnnot",
                                 "pbCHR", "pbPOS", "pbEND", "pbAnnot", 
                                 "ovpLen")
################################################################################
# extract "relevant" pacbio annotations on the overlaps
extractPacbioInfoOnOverlaps <- function(x) {
    info <- strsplit(x$pbAnnot, ";")[[1]]
    ls <- list()
    ls[["pacbioQual"]] <- info[2]
    ls[["pacbioREPEAT"]] <- info[5]
    unlist(ls)
}

temp <- adply(gatkIntersectPacbio, 1, extractPacbioInfoOnOverlaps)
gatkIntersectPacbio <- temp

temp <- adply(mantaIntersectPacbio, 1, extractPacbioInfoOnOverlaps)
mantaIntersectPacbio <- temp

rm(temp, extractPacbioInfoOnOverlaps)
################################################################################
# stratify PacBio records by length
classifyLen <- function(l) {
    if(l<301) {
        1;
    } else if (l<1001) {
        2;
    } else {
        3;
    }
}
factoredPacbioLength <- factor(unlist(lapply(background_pacbio$INFO, 
                                       function(x) classifyLen(as.numeric(sub("SVLEN=", "", strsplit(x, ";")[[1]][3]))))))
rm(classifyLen)
################################################################################
# association with manta QUAL and GQ
# also association on manta's own calls, between SVLEN and QUAL/GQ/GT:
pdf("Results/associationWithPacBioQualAndRepeatTypes.pdf")
par(mfrow=c(2,2))
# Look at types of repeats of the overlapping GATK/PacBio records, Manta/PacBio records,
# and PacBio background records
par(mar=c(5.1,3.0,3.1,0.5))
# par(mar=c(5.1,4.1,4.1,2.1))
# (1,1)
temp <- factor(gatkIntersectPacbio$pacbioREPEAT)
barplot(table(temp), ylim=c(0, 1200),
        las=2, col = rgb(0,0,1,0.5), cex.names = 0.7, 
        main = "\"Validated\" short-read del, by repeat types")
temp <- factor(mantaIntersectPacbio$pacbioREPEAT)
barplot(table(temp), ylim=c(0, 1200),
        las=2, col = rgb(1,0,0,0.5), cex.names = 0.7, 
        main = "", add = T)
legend("topright", 
       legend=c("Manta", "GATK"), 
       fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), 
       text.font=2, cex=0.8)

# (1,2)
temp <- factor(apply(background_pacbio, 1, 
                  function(x)
                      sub("REPEAT_TYPE=", "", str_extract(x[8], "REPEAT_TYPE=[A-Za-z0-9]+(_[A-Za-z0-9]+)?"))))
z <- barplot(table(temp), ylim=c(0, 1200),
             las=2, col = rgb(0,1,0,0.5), xpd = F, cex.names = 0.7, 
             main = "PacBio deletion, by repeat types")
text(x=z[13]-0.5, y = 1150, table(temp)[["NotMasked"]], cex=0.7)
text(x=z[17]-0.5, y = 1150, table(temp)[["TRF"]], cex=0.7)

par(mar=c(5.1,4.1,1.5,2.1))
# (2,1)
qqplot(log10(as.integer(sub("SVLEN=", "", str_extract(background_pacbio$INFO[grep("TRF", background_pacbio$INFO)], "SVLEN=[0-9]+")))), 
       y=c(1:length(background_pacbio$INFO[grep("TRF", background_pacbio$INFO)])), 
       xlab="TRF variant size", 
       ylab="PacBio TRF deletions sorted by size", 
       type='l', xaxt='n')
ticks <- seq(0, ceiling(4), by=1)
ticks2 <- seq(0, ceiling(4), by=0.2)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=ticks2, labels = FALSE)
axis(1, at=ticks, labels = labels)
abline(v = log10(500), lty=2)
abline(v = log10(1000), lty=2)
axis(1, at=log10(500), padj=-1,
     labels = "500",tick = FALSE, cex=0.8)
# (2,2)
qqplot(log10(as.integer(sub("SVLEN=", "", str_extract(background_pacbio$INFO[grep("NotMasked", background_pacbio$INFO)], "SVLEN=[0-9]+")))),
       y=c(1:length(background_pacbio$INFO[grep("NotMasked", background_pacbio$INFO)])), 
       xlab="NotMasked variant size", 
       ylab="PacBio NotMasked deletions sorted by size", 
       type='l', xaxt='n')
ticks <- seq(0, ceiling(5), by=1)
ticks2 <- seq(0, ceiling(5), by=0.2)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=ticks2, labels = FALSE)
axis(1, at=ticks, labels = labels)
abline(v = log10(500), lty=2)
abline(v = log10(1000), lty=2)
axis(1, at=log10(500), padj=-1,
     labels = "500",tick = FALSE, cex=0.8)
dev.off()
# system("pdf2ps associationWithPacBioQualAndRepeatTypes.pdf associationWithPacBioQualAndRepeatTypes.eps")
################################################################################
# for(i in c(1:10)) { ## add HET/HOMVAR numbers
#     qual_low <- histogram$breaks[i]
#     qual_high <- histogram$breaks[i+1]
#     height <- histogram$counts[i]+15
#     a <- length(which(mantaGTs[background_manta$QUAL>qual_low & background_manta$QUAL<=qual_high]=="0/1"))
#     b <- length(which(mantaGTs[background_manta$QUAL>qual_low & background_manta$QUAL<=qual_high]=="1/1"))
#     c <- length(which(overlappedGATKcalls_mantaGTs[overlappedGATKcalls_mantaQuals>qual_low & overlappedGATKcalls_mantaQuals<=qual_high]=="0/1"))
#     d <- length(which(overlappedGATKcalls_mantaGTs[overlappedGATKcalls_mantaQuals>qual_low & overlappedGATKcalls_mantaQuals<=qual_high]=="1/1"))
#     text(paste(a, b, sep="/"), x=qual_high-50, y=height, font=2, cex=0.4)
#     text(paste(c, d, sep="/"), x=qual_high-50, y=height-30, font=2, cex=0.4)
# }
################################################################################
# # Look at the non-linear correlation between the Manata and PacBio overlapping calls.
# cat("", file = T)
# cat("Next look at the Spearman correlation between the", nrow(mantaIntersectPacbio), "overlapping Manta and PacBio calls", sep=" ")
# extractQUALs <- function(x) {
#     ls <- list()
#     ls[["mantaQual"]] <- strsplit(x$mantaAnnot, ";")[[1]][2]
#     ls[["pbQual"]] <- strsplit(x$pbAnnot, ";")[[1]][2]
#     unlist(ls)
# }
# temp <- adply(mantaIntersectPacbio, 1, extractQUALs)
# library(Hmisc)
# rcorr(cbind(as.integer(temp$pbQual), as.integer(temp$mantaQual)), type = "spearman")
# rm(temp, extractQUALs)
################################################################################
# (1,1)
# histogram <- hist(background_pacbio$QUAL, 
#                   breaks = seq(from=0, to=100, by=10),
#                   col=rgb(0,1,0,0.5), 
#                   main="Overlapping and PacBio deletions", 
#                   xlab="QUAL", ylab="counts")
# hist(as.integer(gatkIntersectPacbio$pacbioQual),
#      breaks = seq(from=0, to=100, by=10),
#      col=rgb(0,0,1,0.5),
#      add=T)
# hist(as.integer(mantaIntersectPacbio$pacbioQual),
#      breaks = seq(from=0, to=100, by=10),
#      col=rgb(1,0,0,0.5), 
#      add=T)
# legend(x = 5, y = 4000, 
#        legend=c("PacBio", "GATK", "Manta"),
#        fill = c(rgb(0,1,0,0.5), rgb(0,0,1,0.5), rgb(1,0,0,0.5)),
#        cex=0.7)
# 
# # (1,2)
# hist(background_pacbio$QUAL[factoredPacbioLength==1], 
#      col=rgb(1,0,0,0.5), 
#      breaks = seq(from=0, to=100, by=10),
#      main = "PacBio deletions, stratified by length", 
#      xlab = "QUAL", ylab = "counts")
# hist(background_pacbio$QUAL[factoredPacbioLength==2], 
#      breaks = seq(from=0, to=100, by=10),
#      col=rgb(0,1,0,0.5),
#      add = T)
# hist(background_pacbio$QUAL[factoredPacbioLength==3],
#      breaks = seq(from=0, to=100, by=10),
#      col=rgb(0,0,1,0.5),
#      add = T)
# legend(x = 5, y = 3000,
#        legend=c("small: [50,300]  ", 
#                 "mid:   (300,1000]", 
#                 expression(paste("large: (1000,", infinity, ")"))),
#        fill = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(0,0,1,0.5)), 
#        cex = 0.7)