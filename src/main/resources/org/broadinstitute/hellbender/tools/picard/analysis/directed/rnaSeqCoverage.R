# Script to generate a normalized coverage vs. position along transcript plot.
#
# @author Tim Fennell

# Parse the arguments
args <- commandArgs(trailing = TRUE)
metricsFile  <- args[1]
outputFile   <- args[2]
bamName  <- args[3]
subtitle <- ifelse(length(args) < 4, "", args[4])

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

data <- read.table(metricsFile, header=T, sep="\t", skip=secondBlankLine, check.names=FALSE)

# The histogram has a normalized_position and normalized_coverage column for each metric "level"
# This code parses out the distinct levels so we can output one graph per level
headers <- sapply(sub(".normalized_coverage","",names(data),fixed=TRUE), "[[" ,1)

## Duplicated header names cause this to barf. KT & Yossi report that this is going to be extremely difficult to
## resolve and it's unlikely that anyone cares anyways. Trap this situation and avoid the PDF so it won't cause
## the workflow to fail
if (any(duplicated(headers))) {
  print(paste("Not creating insert size PDF as there are duplicated header names:", headers[which(duplicated(headers))]))
} else {
    pdf(outputFile)
    levels <- c()
    for (i in 2:length(headers)) {
        if (!(headers[i] %in% levels)) {
            levels[length(levels)+1] <- headers[i]
        }
    }

    # Some constants that are used below
    COLORS = c("royalblue", "#FFAAAA", "palegreen3");

    # For each level, plot of the normalized coverage by GC
    for (i in 1:length(levels)) {

        # Reconstitutes the histogram column header for this level
        nc <- paste(levels[i], "normalized_coverage", sep=".")

        plot(x=data$normalized_position, y=as.matrix(data[nc]),
            type="o",
            xlab="Normalized Distance Along Transcript",
            ylab="Normalized Coverage",
            xlim=c(0, 100),
            ylim=range(0, max(data[nc])),
            col="royalblue",
            main=paste("RNA-Seq Coverage vs. Transcript Position\n", levels[i], " ", ifelse(subtitle=="", "", paste("(", subtitle, ")", sep="")), "\nin file ", bamName,sep=""))

        # Add a horizontal line at coverage=1
        abline(h=1, col="lightgrey");
    }
    dev.off();
}