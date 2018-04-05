# Collects histogram, QQ-plot (against uniform distribution) of insertion 
# and deletion calls

#########################################################################

args<-commandArgs(TRUE)
wd <- args[1]
wd <- paste0(wd, "/Accounting")
setwd(wd)
print(paste("working in directory", wd))

mantaAndPacBioDir <- args[2]

scriptDir <- args[3]
helperPath <- paste0(scriptDir, "AccountingAndPrep/func_collectSizes.R")
source( helperPath )

sizesCollection <- collectSizes(paste0(mantaAndPacBioDir, "/mantaExactCleanDeletionSizes.txt"), 
                                paste0(mantaAndPacBioDir, "/mantaExactCleanInsertionSizes.txt"),
                                paste0(mantaAndPacBioDir, "/mantaExactCleanTandupSizes.txt"),
                                "gatkCleanDeletionSizes.txt",
                                "gatkInsertionSizes.txt",
                                "gatkTandupSizes.txt",
                                paste0(mantaAndPacBioDir, "/pacbioDeletionSizes_0.txt"),
                                paste0(mantaAndPacBioDir, "/pacbioInsertionSizes_0.txt"), 
                                paste0(mantaAndPacBioDir, "/InsDupRPL/Manta_PASS_PRECISE_nonBND_primaryContigs_longRangeSubstitution.txt"),
                                "../InsDupRPL/GATK_primaryContigs_scarredDel.txt")
#########################################################################
pdf("insDelHistogramsAndQQplot.pdf", onefile = F)
layout(matrix(c(1,1,1,2,3,4,5,6,7), 3,3, byrow = TRUE), 
       heights = c(1.5,1,1,1,1,1,1))
colors <- c("black", "magenta3", "red", "orange", "blue", "green")
x_limit <- log10(max(c(sizesCollection$cleanDelSizes_manta, 
                       sizesCollection$cleanDelSizes_gatk, 
                       sizesCollection$insertionAndDup_manta, 
                       sizesCollection$insertionAndDup_gatk, 
                       sizesCollection$pacbioDeletionSizes, 
                       sizesCollection$pacbioInsertionSizes)))
y_limit <- max(c(length(sizesCollection$cleanDelSizes_manta), 
                 length(sizesCollection$cleanDelSizes_gatk), 
                 length(sizesCollection$insertionAndDup_manta), 
                 length(sizesCollection$insertionAndDup_gatk), 
                 length(sizesCollection$pacbioDeletionSizes), 
                 length(sizesCollection$pacbioInsertionSizes)))

########### QQ plot of log(variant sizes) against uniform distribution 
# row 1
qqplot(x=log10(sizesCollection$pacbioInsertionSizes), 
       y=1:length(sizesCollection$pacbioInsertionSizes), 
       xlim=c(log10(49), x_limit), 
       ylim=c(0, y_limit),
       xlab="variant size", 
       ylab="variant discovered sorted by size", 
       main = "Q-Q plot of variants size vs variants discovered",
       xaxt='n',type='l',col=colors[1])
lines(x=log10(sizesCollection$pacbioDeletionSizes), 
      y=1:length(sizesCollection$pacbioDeletionSizes), 
      type='l', col=colors[2])

lines(x=log10(sizesCollection$cleanDelSizes_manta), 
      y=1:length(sizesCollection$cleanDelSizes_manta), 
      type='l', col=colors[3])
lines(x=log10(sizesCollection$insertionAndDup_manta), 
      y=1:length(sizesCollection$insertionAndDup_manta), 
      type='l', col=colors[4])
lines(x=log10(sizesCollection$cleanDelSizes_gatk), 
      y=1:length(sizesCollection$cleanDelSizes_gatk), 
      type='l', col=colors[5])
lines(x=log10(sizesCollection$insertionAndDup_gatk), 
      y=1:length(sizesCollection$insertionAndDup_gatk), 
      type='l', col=colors[6])

legend("topright", 
       legend=c("PacBio ins", "PacBio del", 
                "Manta clean del", "Manta ins+dup", 
                "GATK clean del", "GATK ins+dup"), 
       col=colors, lty = 1, text.font=2, lwd=3)

ticks <- seq(0, ceiling(x_limit), by=1)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=ticks, labels=labels)

# some decoration
abline(v=log10(250), lty=2)
abline(v=log10(350), lty=3)
abline(v=log10(5500), lty=2)
abline(v=log10(6500), lty=3)

axis(1, at=c(log10(250), log10(5500)), 
     labels = c("250", "5500"),
     padj=-1.5, tick = FALSE, cex=0.8)
axis(3, at=c(log10(350), log10(6500)),
     labels = c("350", "6500"), 
     padj=1.5, tick = FALSE, cex=0.8)


########### historgram
par(mar=c(3.1,3.0,1.1,0.5))

## row 2: for SINE
# deletion
hist(sizesCollection$pacbioDeletionSizes[sizesCollection$pacbioDeletionSizes<1001], 
     breaks = seq(from=50, to=1000, by=50), 
     col=rgb(0,1,0,0.5), main="clean deletions", xlab="", ylab="")
hist(sizesCollection$cleanDelSizes_manta[sizesCollection$cleanDelSizes_manta<1001], 
     breaks = seq(from=50, to=1000, by=50), 
     col=rgb(1,0,0,0.5), add=T)
hist(sizesCollection$cleanDelSizes_gatk[sizesCollection$cleanDelSizes_gatk<1001], 
     breaks = c(seq(from=50, to=1000, by=50)), 
     col=rgb(0,0,1,0.5), add=T)
legend("topright", 
       legend=c("PacBio", "Manta", "GATK"), 
       fill = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), rgb(0,0,1,0.5)),
       text.font=2, cex=0.7)
# insertion
hist(sizesCollection$pacbioInsertionSizes[sizesCollection$pacbioInsertionSizes<1001], 
     breaks = seq(from=50, to=1000, by=50), 
     col=rgb(0,1,0,0.5), main="ins + tanDup", xlab="", ylab="")
hist(sizesCollection$insertionAndDup_manta[sizesCollection$insertionAndDup_manta<1001], 
     breaks = seq(from=50, to=1000, by=50), 
     col=rgb(1,0,0,0.5), add=T)
hist(sizesCollection$insertionAndDup_gatk[sizesCollection$insertionAndDup_gatk<1001], 
     breaks = c(seq(from=50, to=1000, by=50)),
     col=rgb(0,0,1,0.5), add=T)
# scared deletion
hist(-sizesCollection$scarDel_manta_del[sizesCollection$scarDel_manta_del<=500], 
     xlim = c(-500,500),
     breaks = seq(from=-500, to=500, by=50), 
     col=rgb(1,0,0,0.5), 
     main="scarred deletions", xlab="", ylab="")
hist(sizesCollection$scarDel_manta_ins[sizesCollection$scarDel_manta_ins<=500], 
     xlim = c(-500,500),
     breaks = seq(from=-500, to=500, by=50), 
     col=rgb(1,0,0,0.5), 
     add = T)
hist(-sizesCollection$scarDel_gatk_del[sizesCollection$scarDel_gatk_del<=500], 
     xlim = c(-500,500),
     breaks = seq(from=-500, to=500, by=50), 
     col=rgb(0,0,1,0.5), add = T)
hist(sizesCollection$scarDel_gatk_ins[sizesCollection$scarDel_gatk_ins<=500], 
     xlim = c(-500,500),
     breaks = seq(from=-500, to=500, by=50), 
     col=rgb(0,0,1,0.5), 
     add = T)

## row 3: for LINE
# deletion
hist(sizesCollection$pacbioDeletionSizes[sizesCollection$pacbioDeletionSizes>1000 & sizesCollection$pacbioDeletionSizes<10000], 
     breaks = seq(from=1000, to=10000, by=500),
     col=rgb(0,1,0,0.5), main="", xlab="", ylab="")
hist(sizesCollection$cleanDelSizes_manta[sizesCollection$cleanDelSizes_manta>1000 & sizesCollection$cleanDelSizes_manta<10000], 
     breaks = seq(from=1000, to=10000, by=500),
     col=rgb(1,0,0,0.5), add=T)
hist(sizesCollection$cleanDelSizes_gatk[sizesCollection$cleanDelSizes_gatk>1000 & sizesCollection$cleanDelSizes_gatk<10000], 
     breaks = seq(from=1000, to=10000, by=500),
     col=rgb(0,0,1,0.5), add=T)
# insertion
hist(sizesCollection$pacbioInsertionSizes[sizesCollection$pacbioInsertionSizes>1000 & sizesCollection$pacbioInsertionSizes<10000], 
     breaks = seq(from=1000, to=10000, by=500),
     col=rgb(0,1,0,0.5), main="", xlab="", ylab="")
hist(sizesCollection$insertionAndDup_manta[sizesCollection$insertionAndDup_manta>1000 & sizesCollection$insertionAndDup_manta<10000], 
     breaks = seq(from=1000, to=10000, by=500),
     col=rgb(1,0,0,0.5), add=T)
hist(sizesCollection$insertionAndDup_gatk[sizesCollection$insertionAndDup_gatk>1000 & sizesCollection$insertionAndDup_gatk<10000], 
     breaks = seq(from=1000, to=10000, by=500),
     col=rgb(0,0,1,0.5), add=T)
# scarred deletion
hist(-sizesCollection$scarDel_manta_del[sizesCollection$scarDel_manta_del>1000 &
                                            sizesCollection$scarDel_manta_del<10000], 
     xlim = c(-10000,10000),
     ylim = c(0,50),
     breaks = seq(from=-10000, to=10000, by=1000), 
     col=rgb(1,0,0,0.5), 
     main="", xlab="", ylab="")
hist(sizesCollection$scarDel_manta_ins[sizesCollection$scarDel_manta_ins>1000 &
                                           sizesCollection$scarDel_manta_ins<10000], 
     xlim = c(-1000,1000),
     ylim = c(0,50),
     breaks = seq(from=-10000, to=10000, by=1000), 
     col=rgb(1,0,0,0.5), 
     add = T)
hist(-sizesCollection$scarDel_gatk_del[sizesCollection$scarDel_gatk_del>1000 &
                                           sizesCollection$scarDel_gatk_del<10000], 
     xlim = c(-1000,1000),
     ylim = c(0,50),
     breaks = seq(from=-10000, to=10000, by=1000), 
     col=rgb(0,0,1,0.5), add = T)
hist(sizesCollection$scarDel_gatk_ins[sizesCollection$scarDel_gatk_ins>1000 &
                                          sizesCollection$scarDel_gatk_ins<1000], 
     xlim = c(-1000,1000),
     ylim = c(0,50),
     breaks = seq(from=-10000, to=10000, by=1000), 
     col=rgb(0,0,1,0.5), 
     add = T)
########### finish up
dev.off()
# system("pdf2ps insDelHistogramsAndQQplot.pdf insDelHistogramsAndQQplot.eps")
