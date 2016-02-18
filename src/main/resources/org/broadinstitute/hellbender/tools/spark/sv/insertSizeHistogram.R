## script to generate histogram of insert sizes from metrics file
## expecting 3 arguments:
## first is the metrics file with the histogram info
## second is the output file
## third is a name for the plot

args <- commandArgs(trailing=TRUE)
metricsFile <- args[1]
pdfFile <- args[2]
bamName <- args[3]
histoWidth <- ifelse(length(args) < 4, 0, as.numeric(args[4]))

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

histogram <- read.table(metricsFile, header=TRUE, sep="\t", skip=secondBlankLine, comment.char="", quote='', check.names=FALSE)

## The histogram has a fr_count/rf_count/tandem_count for each metric "level"
## This code parses out the distinct levels so we can output one graph per level
headers <- sapply(sub(".fr_count","",names(histogram),fixed=TRUE), "[[" ,1)
headers <- sapply(sub(".rf_count","",headers,fixed=TRUE), "[[" ,1)
headers <- sapply(sub(".tandem_count","",headers,fixed=TRUE), "[[" ,1)

## Duplicated header names cause this to barf. KT & Yossi report that this is going to be extremely difficult to
## resolve and it's unlikely that anyone cares anyways. Trap this situation and avoid the PDF so it won't cause
## the workflow to fail
if (any(duplicated(headers))) {
  print(paste("Not creating insert size PDF as there are duplicated header names:", headers[which(duplicated(headers))]))
} else {
  levels <- c()
  for (i in 2:length(headers)) {
    if (!(headers[i] %in% levels)) {
      levels[length(levels)+1] <- headers[i]
    }
  }

  pdf(pdfFile)

  for (i in 1:length(levels)) {
    ## Reconstitutes the histogram column headers for this level
    fr <- paste(levels[i], "fr_count", sep=".")
    rf <- paste(levels[i], "rf_count", sep=".")
    tandem <- paste(levels[i], "tandem_count", sep=".")

    frrange = ifelse(fr %in% names(histogram), max(histogram[fr]), 0)
    rfrange = ifelse(rf %in% names(histogram), max(histogram[rf]), 0)
    tandemrange = ifelse(tandem %in% names(histogram), max(histogram[tandem]), 0)

    yrange <- max(frrange, rfrange, tandemrange)
    xrange <- ifelse(histoWidth > 0, histoWidth, max(histogram$insert_size))

    plot(x=NULL, y=NULL,
         type="n",
         main=paste("Insert Size Histogram for", levels[i], "\nin file", bamName),
         xlab="Insert Size",
         ylab="Count",
         xlim=range(0, xrange),
         ylim=range(0, yrange))

    colors <- c()
    labels <- c()

    if (fr %in% names(histogram) ) {
      lines(histogram$insert_size, as.matrix(histogram[fr]),  type="h", col="red")
      colors <- c(colors, "red")
      labels <- c(labels, "FR")
    }
    if (rf %in% names(histogram)) {
      lines(histogram$insert_size, as.matrix(histogram[rf]),  type="h", col="blue")
      colors <- c(colors, "blue")
      labels <- c(labels, "RF")
    }

    if (tandem %in% names(histogram)) {
      lines(histogram$insert_size, as.matrix(histogram[tandem]),  type="h", col="orange")
      colors <- c(colors, "orange")
      labels <- c(labels, "TANDEM")
    }

    ## Create the legend
    legend("topright", labels, fill=colors, col=colors, cex=0.7)
  }

  dev.off()
}

