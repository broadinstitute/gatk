# Script to generate a chart of quality score distribution in a file
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

plot(histogram$QUALITY,
     histogram$COUNT_OF_Q,
     type="n",
     main=paste("Quality Score Distribution\nin file ",bamFile," ",ifelse(subtitle == "","",paste("(",subtitle,")",sep="")),sep=""),
     xlab="Quality Score",
     ylab="Observations")

qColor  <- "blue"
oqColor <- "lightcyan2"
width <- 5

# Plot OQ first so that it's "behind" the regular qualities
if (!is.null(histogram$COUNT_OF_OQ)) {
    lines(histogram$QUALITY+0.25, histogram$COUNT_OF_OQ, type="h", col=oqColor, lty=1, lwd=width, lend="square");
}

# Then plot the regular qualities
lines(histogram$QUALITY, histogram$COUNT_OF_Q, type="h", col=qColor, lty=1, lwd=width, lend="square");

# And add a legend
legend("topleft", pch=c(15,15), legend=c("Quality Scores", "Original Quality Scores"), col=c(qColor, oqColor))

dev.off()

