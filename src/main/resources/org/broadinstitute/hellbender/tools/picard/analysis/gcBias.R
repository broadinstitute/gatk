# Script to generate a chart to display GC bias based upon read starts observed
# in windows along the genome.
#
# @author Tim Fennell

# Parse the arguments
args <- commandArgs(trailing=T)
metricsFile  <- args[1]
outputFile   <- args[2]
datasetName  <- args[3]
subtitle  <- args[4]
windowSize   <- args[5]

# Figure out where the metrics and the histogram are in the file and parse them out
startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)

firstBlankLine=0

for (i in 1:length(startFinder)) {
        if (startFinder[i] == "") {
                if (firstBlankLine==0) {
                        firstBlankLine=i+1
                } else {
                        secondBlankLine=i+1
                        break
                }
        }
}

metrics <- read.table(metricsFile, header=T, sep="\t", skip=firstBlankLine)
pdf(outputFile)

# Some constants that are used below
Y_AXIS_LIM = 2;
MAX_QUALITY_SCORE = 40;
COLORS = c("royalblue", "#FFAAAA", "palegreen3");

# Adjust to give more margin on the right hand side
par(mar = c(5, 4, 4, 4));

# Do the main plot of the normalized coverage by GC
plot(type="p", x=metrics$GC, y=metrics$NORMALIZED_COVERAGE,
     xlab=paste(c("GC% of", windowSize, "base windows"), sep=" ", collapse=" "),
     ylab="Fraction of normalized coverage",
     xlim=c(0,100),
     ylim=c(0, Y_AXIS_LIM),
     col=COLORS[1],
     main=paste(datasetName, "GC Bias Plot", "\n", subtitle)
    );

# Add lines at the 50% GC and coverage=1
abline(h=1, v=50, col="lightgrey");

# Add error bars
arrows(metrics$GC,
       metrics$NORMALIZED_COVERAGE - metrics$ERROR_BAR_WIDTH,
       metrics$GC,
       metrics$NORMALIZED_COVERAGE + metrics$ERROR_BAR_WIDTH,
       code = 3, angle = 90, length = 0.05, col="grey");

# Plot count of windows as a separate series near the bottom
window_ratio = 0.5 / max(metrics$WINDOWS);
scaled_windows = metrics$WINDOWS * window_ratio;
lines(metrics$GC, scaled_windows, type="h", col=COLORS[2], lwd=3);

# Plot the quality series
lines(metrics$GC, metrics$MEAN_BASE_QUALITY * Y_AXIS_LIM / MAX_QUALITY_SCORE, type="l", col=COLORS[3]);
axis(4,
     at=c(0, Y_AXIS_LIM/4, Y_AXIS_LIM/4*2, Y_AXIS_LIM/4*3, Y_AXIS_LIM),
     labels=c(0, MAX_QUALITY_SCORE/4, MAX_QUALITY_SCORE/4*2, MAX_QUALITY_SCORE/4*3, MAX_QUALITY_SCORE)
    );
mtext("Mean base quality", side=4, line=2.5);

# And finally add a legend
legend("topleft", pch=c(1,15, 45), legend=c("Normalized Coverage", "Windows at GC%", "Base Quality at GC%"), col=COLORS)

dev.off();