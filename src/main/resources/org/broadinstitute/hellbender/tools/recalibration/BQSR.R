library("ggplot2")
library(gplots)
library("reshape")
library("grid")
library("tools") #For compactPDF in R 2.13+
library(gsalib)



if ( interactive() ) {
  args <- c("NA12878.6.1.dedup.realign.recal.bqsr.grp.csv", "NA12878.6.1.dedup.realign.recal.bqsr.grp", NA)
} else {
  args <- commandArgs(TRUE)
} 

data <- read.csv(args[1])

data$Recalibration = as.factor(sapply(as.character(data$Recalibration),function(x) { 
  xu = toupper(x); 
  if (xu == "ORIGINAL") "BEFORE" else 
    if (xu == "RECALIBRATED") "AFTER" else 
      if (xu == "RECALIBRATION") "BQSR" else
        xu }));

gsa.report <- gsa.read.gatkreport(args[2])

gsa.report$Arguments$Value = as.character(gsa.report$Arguments$Value);
gsa.report$Arguments = subset(gsa.report$Arguments,subset= Argument != "plot_pdf_file");
if (length(levels(data$Recalibration)) > 1) {
  gsa.report$Arguments = subset(gsa.report$Arguments,subset= Argument != "recalibration_report");
}
gsa.report$Arguments$Value[gsa.report$Argument$Value == "null"] = "None";

gsa.report.covariate.argnum = gsa.report$Arguments$Argument == "covariate";
gsa.report$Arguments$Value[gsa.report.covariate.argnum] = sapply(strsplit(gsa.report$Arguments$Value[gsa.report.covariate.argnum],","),function(x) { 
  y = sub("(^.+)Covariate","\\1",x); paste(y,collapse=",") } );

data <- within(data, EventType <- factor(EventType, levels = rev(levels(EventType))))

numRG = length(unique(data$ReadGroup))
blankTheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks = element_blank())

# Viewport (layout 2 graphs top to bottom)
distributeGraphRows <- function(graphs, heights = c()) {
  if (length(heights) == 0) {
    heights <- rep.int(1, length(graphs))
  }
  heights <- heights[!is.na(graphs)]
  graphs <- graphs[!is.na(graphs)]
  numGraphs <- length(graphs)
  Layout <- grid.layout(nrow = numGraphs, ncol = 1, heights=heights)
  grid.newpage()
  pushViewport(viewport(layout = Layout))
  subplot <- function(x) viewport(layout.pos.row = x, layout.pos.col = 1)
  for (i in 1:numGraphs) {
    print(graphs[[i]], vp = subplot(i))
  }
}


for(cov in levels(data$CovariateName)) {    # for each covariate in turn  
  d = data[data$CovariateName==cov,]        # pull out just the data for this covariate so we can treat the non-numeric values appropriately
  if( cov == "Context" ) {
    d$CovariateValue = as.character(d$CovariateValue)
    d$CovariateValue = substring(d$CovariateValue,nchar(d$CovariateValue)-2,nchar(d$CovariateValue))
  } else {
    d$CovariateValue = as.numeric(levels(d$CovariateValue))[as.integer(d$CovariateValue)] # efficient way to convert factors back to their real values
  }
  #d=subset(d,Observations>2000) # only show bins which have enough data to actually estimate the quality
  dSub=subset(d,EventType=="Base Substitution")
  dIns=subset(d,EventType=="Base Insertion")
  dDel=subset(d,EventType=="Base Deletion")
  dSub=dSub[sample.int(length(dSub[,1]),min(length(dSub[,1]),2000)),] # don't plot too many values because it makes the PDFs too massive
  dIns=dIns[sample.int(length(dIns[,1]),min(length(dIns[,1]),2000)),] # don't plot too many values because it makes the PDFs too massive
  dDel=dDel[sample.int(length(dDel[,1]),min(length(dDel[,1]),2000)),] # don't plot too many values because it makes the PDFs too massive
  d=rbind(dSub, dIns, dDel)
  
  if( cov != "QualityScore" ) {    
    p <- ggplot(d, aes(x=CovariateValue,y=Accuracy,alpha=log10(Observations))) + ylim(min(-10,d$Accuracy),max(10,d$Accuracy)) +
      geom_abline(intercept=0, slope=0, linetype=2) + 
      xlab(paste(cov,"Covariate")) +
      ylab("Quality Score Accuracy") +
      blankTheme
    if(cov == "Cycle") {
      b <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("BEFORE"="maroon1","AFTER"="blue","BQSR"="black")) + facet_grid(.~EventType) +
        theme(axis.text.x=element_text(angle=90, hjust=0))
      
      p <- ggplot(d, aes(x=CovariateValue,y=AverageReportedQuality,alpha=log10(Observations))) +
        xlab(paste(cov,"Covariate")) +
        ylab("Mean Quality Score") + ylim(0,max(42,d$AverageReportedQuality)) +
        blankTheme
      e <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("BEFORE"="maroon1","AFTER"="blue","BQSR"="black")) + facet_grid(.~EventType) +
        theme(axis.text.x=element_text(angle=90, hjust=0))
      
      
    } else {
      c <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("BEFORE"="maroon1","AFTER"="blue","BQSR"="black")) + facet_grid(.~EventType) +
        theme(axis.text.x=element_text(angle=90, hjust=0)) + xlab(paste(cov,"Covariate (3 base suffix)"))
      p <- ggplot(d, aes(x=CovariateValue,y=AverageReportedQuality,alpha=log10(Observations))) +
        xlab(paste(cov,"Covariate (3 base suffix)")) +
        ylab("Mean Quality Score") +
        blankTheme
      f <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("BEFORE"="maroon1","AFTER"="blue","BQSR"="black")) + facet_grid(.~EventType) +
        theme(axis.text.x=element_text(angle=90, hjust=0))
      
    }
  } else {
    p <- ggplot(d, aes(x=AverageReportedQuality,y=EmpiricalQuality,alpha=log10(Observations))) +
      geom_abline(intercept=0, slope=1, linetype=2) + 
      xlab("Reported Quality Score") +
      ylab("Empirical Quality Score") +
      blankTheme
    a <- p + geom_point(aes(color=Recalibration)) + scale_color_manual(values=c("BEFORE"="maroon1","AFTER"="blue","BQSR"="black")) + facet_grid(.~EventType)
    
    p <- ggplot(d, aes(x=CovariateValue)) +
      xlab(paste(cov,"Covariate")) +
      ylab("No. of Observations (area normalized)") +
      blankTheme
    d <- p + geom_histogram(aes(fill=Recalibration,weight=Observations,y=..ndensity..),alpha=0.6,binwidth=1,position="identity")
    d <- d + scale_fill_manual(values=c("BEFORE"="maroon1","AFTER"="blue","BQSR"="black"))
    d <- d + facet_grid(.~EventType) 
    #    d <- d + scale_y_continuous(formatter="comma")
  }
}

if ( ! is.na(args[3]) )
  pdf(args[3],height=9,width=15)

#frame()
textplot(gsa.report$Arguments, show.rownames=F)
title(
  main="GATK BaseRecalibration report",
  sub=date())

distributeGraphRows(list(a,b,c), c(1,1,1))
distributeGraphRows(list(d,e,f), c(1,1,1))

# format the overall information
rt0 <- data.frame(
  ReadGroup = gsa.report$RecalTable0$ReadGroup,
  EventType = gsa.report$RecalTable0$EventType,
  EmpiricalQuality = sprintf("%.1f", gsa.report$RecalTable0$EmpiricalQuality),
  EstimatedQReported = sprintf("%.1f", gsa.report$RecalTable0$EstimatedQReported),
  Observations = sprintf("%.2e", gsa.report$RecalTable0$Observations),
  Errors = sprintf("%.2e", gsa.report$RecalTable0$Errors))  
textplot(t(rt0), show.colnames=F)
title("Overall error rates by event type")

# plot per quality score recalibration table
textplot(gsa.report$RecalTable1, show.rownames=F)
title("Error rates by event type and initial quality score")

if ( ! is.na(args[3]) ) {
  dev.off()
  if (exists('compactPDF')) {
    compactPDF(args[2])
  }
}

