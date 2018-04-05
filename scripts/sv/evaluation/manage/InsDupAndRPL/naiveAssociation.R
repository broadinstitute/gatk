## This script tries to look at, naively, what type of variables are associated 
##  with the overlapping calls between PacBio and GATK, and what's associated with 
##  those that are present in PacBio's call set but not in GATK's call set, 
##  and vice versa.
library("plotrix")
library("plyr")
library("stringr")
################################################################################
# load two background data, manta call set and gatk call set
#   and load intersections
background_pacbio <- read.table("PacBio_primaryContigs_ins.vcf", 
                                stringsAsFactors = F, header = F)
names(background_pacbio) <- c("CHROM", "POS", "ID", "REF", "ALT", 
                              "QUAL", "FILTER", "INFO")
################################################################################
# intersections
gatkIntersectPacBio_insIns <- read.table("Results/insIns.GATKvsPacBio.window_w50.txt", 
                                         stringsAsFactors = F, header = F, 
                                         na.strings = c(".", "-1", "0"))
gatkIntersectPacBio_dupIns <- read.table("Results/dupIns.GATKvsPacBio.window_l0r50.txt", 
                                         stringsAsFactors = F, header = F, 
                                         na.strings = c(".", "-1", "0"))
names(gatkIntersectPacBio_insIns) <- c("gatkCHR", "gatkPOS", "gatkID", "gatkREF", "gatkALT", "gatkQUAL", 
                                       "gatkFILTER", "gatkINFO", 
                                       "pbCHR", "pbPOS", "pbID", "pbREF", "pbALT", "pbQUAL", 
                                       "pbFILTER", "pbINFO")
names(gatkIntersectPacBio_dupIns) <- names(gatkIntersectPacBio_insIns)
gatkIntersectPacBio <- rbind(gatkIntersectPacBio_insIns, gatkIntersectPacBio_dupIns)
rm(gatkIntersectPacBio_insIns, gatkIntersectPacBio_dupIns)

mantaIntersectPacBio_insIns <- read.table("Results/insIns.MantavsPacBio.window_w50.txt", 
                                         stringsAsFactors = F, header = F, 
                                         na.strings = c(".", "-1", "0"))
mantaIntersectPacBio_dupIns <- read.table("Results/dupIns.MantavsPacBio.window_l0r50.txt", 
                                         stringsAsFactors = F, header = F, 
                                         na.strings = c(".", "-1", "0"))
names(mantaIntersectPacBio_insIns) <- c("mantaCHR", "mantaPOS", "mantaID", "mantaREF", "mantaALT", "mantaQUAL", 
                                       "mantaFILTER", "mantaINFO", "mantaFORMAT", "mantaSAMPLE",
                                       "pbCHR", "pbPOS", "pbID", "pbREF", "pbALT", "pbQUAL", 
                                       "pbFILTER", "pbINFO")
names(mantaIntersectPacBio_dupIns) <- names(mantaIntersectPacBio_insIns)
mantaIntersectPacBio <- rbind(mantaIntersectPacBio_insIns, mantaIntersectPacBio_dupIns)
rm(mantaIntersectPacBio_insIns,mantaIntersectPacBio_dupIns)
################################################################################
# extract "relevant" pacbio annotations on the overlaps
extractPacBioInfoOnOverlaps <- function(x) {
    ls <- list()
    str_extract(x$pbINFO, "REPEAT_TYPE=[A-Za-z0-9]+(_[A-Za-z0-9]+)?")
    ls[["pacbioREPEAT"]] <- sub("REPEAT_TYPE=", "", str_extract(x$pbINFO, "REPEAT_TYPE=[A-Za-z0-9]+(_[A-Za-z0-9]+)?"))
    unlist(ls)
}

temp <- adply(gatkIntersectPacBio, 1, extractPacBioInfoOnOverlaps)
gatkIntersectPacBio <- temp

temp <- adply(mantaIntersectPacBio, 1, extractPacBioInfoOnOverlaps)
mantaIntersectPacBio <- temp

rm(temp, extractPacBioInfoOnOverlaps)
################################################################################
# association with PacBio repeat type
pdf("Results/associationWithPacBioQualAndRepeatTypes.pdf")
par(mfrow=c(2,2))

pacbioRepeatTypes <- factor(apply(background_pacbio, 1, 
                                  function(x)
                                      sub("REPEAT_TYPE=", "", str_extract(x[8], "REPEAT_TYPE=[A-Za-z0-9]+(_[A-Za-z0-9]+)?"))))

# Look at types of repeats of the overlapping GATK/PacBio records, Manta/PacBio records,
# and PacBio background records
par(mar=c(5.1,3.0,3.1,0.5))
# (1,1)
temp <- table(factor(gatkIntersectPacBio$pacbioREPEAT))
missingRepeats <- setdiff(names(table(pacbioRepeatTypes)), names(temp))
height <- as.vector(temp)
names(height) <- names(temp)
for(name in missingRepeats) {
    height[[name]] <- 0
}
height <- height[order(names(height))]
barplot(height, ylim=c(0, 1200), names.arg = names(height),
        las=2, col = rgb(0,0,1,0.5), cex.names = 0.7, 
        main = "\"Validated\" short-read ins, by repeat types")

temp <- table(factor(mantaIntersectPacBio$pacbioREPEAT))
missingRepeats <- setdiff(names(table(pacbioRepeatTypes)), names(temp))
height <- as.vector(temp)
names(height) <- names(temp)
for(name in missingRepeats) {
    height[[name]] <- 0
}
height <- height[order(names(height))]
barplot(height, ylim=c(0, 1200),
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
             main = "PacBio ins., by repeat types")
text(x=z[13]-0.5, y = 1150, table(temp)[["NotMasked"]], cex=0.7)
text(x=z[17]-0.5, y = 1150, table(temp)[["TRF"]], cex=0.7)

par(mar=c(5.1,4.1,1.5,2.1))
# (2,1)
qqplot(log10(as.integer(sub("SVLEN=", "", str_extract(background_pacbio$INFO[grep("TRF", background_pacbio$INFO)], "SVLEN=[0-9]+")))), 
       y=c(1:length(background_pacbio$INFO[grep("TRF", background_pacbio$INFO)])), 
       xlab="TRF variant size", 
       ylab="PacBio TRF ins. sorted by size", 
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
       ylab="PacBio NotMasked ins. sorted by size", 
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