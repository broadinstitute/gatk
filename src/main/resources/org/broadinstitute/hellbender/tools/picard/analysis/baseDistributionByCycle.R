# Script to generate a chart of the base distribution by cycle
# @author Nils Homer

# Parse the arguments
args <- commandArgs(trailing=T);
metricsFile  <- args[1];
outputFile   <- args[2];
bamFile  <- args[3];
subtitle <- ifelse(length(args) < 4, "", args[4]);


# Figure out where the metrics and the histogram are in the file and parse them out
startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE);

firstBlankLine=0;

for (i in 1:length(startFinder)) {
        if (startFinder[i] == "") {
                if (firstBlankLine==0) {
                        firstBlankLine=i+1;
                } else {
                        secondBlankLine=i+1;
                        break;
                }
        }
}

metrics <- read.table(metricsFile, header=T, sep="\t", skip=firstBlankLine);

# Then plot the histogram as a PDF
pdf(outputFile);

plot(x=c(1, 20+nrow(metrics)),
     y=c(0, max(metrics[,3:7])),
     main=paste("Base Distribution by Cycle\nin file ",bamFile," ",ifelse(subtitle == "","",paste("(",subtitle,")",sep="")),sep=""),
     xlab="Cycle",
     ylab="Base Percentage",
     type="n");

colors = c("red", "orange", "blue", "purple", "black");

for (i in 1:5) {
    lines(x=1:nrow(metrics),
        y=metrics[,2+i],
        col=colors[i],
        type="l",
        lty=1);
}

legend("bottomright", lwd=1, legend=c("PCT_A", "PCT_C", "PCT_G", "PCT_T", "PCT_N"), col=colors);

dev.off();
