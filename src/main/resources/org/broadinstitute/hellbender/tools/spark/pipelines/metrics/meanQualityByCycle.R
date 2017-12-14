# Script to generate a chart of mean quality by cycle from a BAM file
# @author Tim Fennell

# Parse the arguments
args <- commandArgs(trailing=T)
metricsFile  <- args[1]
outputFile   <- args[2]
bamFile  <- args[3]
subtitle <- ifelse(length(args) < 4, "", args[4])


# Figure out where the metrics and the histogram are in the file and parse them out
startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)

firstBlankLine=0

for (i in 1:length(startFinder))
{
        if (startFinder[i] == "") {
                if (firstBlankLine==0) {
                        firstBlankLine=i+1
                } else {
                        secondBlankLine=i+1
                        break
                }
        }
}

metrics <- read.table(metricsFile, header=T, nrows=1, sep="\t", skip=firstBlankLine)
histogram <- read.table(metricsFile, header=T, sep="\t", skip=secondBlankLine)

# Then plot the histogram as a PDF
pdf(outputFile)

plot(histogram$CYCLE,
     histogram$MEAN_QUALITY,
     type="n",
     main=paste("Quality by Cycle\nin file ",bamFile," ",ifelse(subtitle == "","",paste("(",subtitle,")",sep="")),sep=""),
     xlab="Cycle",
     ylab="Mean Quality",
     ylim=range(0,50))

qColor  <- "darkblue"
oqColor <- rgb(1, 0.25, 0.25, 0.75)

# Plot OQ first so that it's "behind" the regular qualities
if (!is.null(histogram$MEAN_ORIGINAL_QUALITY)) {
    lines(histogram$CYCLE, histogram$MEAN_ORIGINAL_QUALITY, type="l", col=oqColor, lty=1);
}

# Then plot the regular qualities
lines(histogram$CYCLE, histogram$MEAN_QUALITY, type="h", col=qColor, lty=1);

# And add a legend
legend("topleft", pch=c(15,15), legend=c("Mean Quality", "Mean Original Quality"), col=c(qColor, oqColor))


dev.off()

