library(dplyr)
library(ggplot2)
library(reshape2)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

round_digits <- -2
files <- list.files(pattern = "summary\\.csv$")
dlist <- lapply(files, read.csv)
names <- lapply(files, function(x) gsub("happy_", "", gsub(".summary.csv", "", x)))
dnamed <- mapply(cbind, dlist, "Name"=names, SIMPLIFY=F)
merged <- Reduce(function(...) merge(..., all=T), dnamed)

names(merged) <- c( "Type", "Filter", "Total", "True Positives", "False Negatives", "QTotal", "False Positives", "Unknown", "Genotype Error", "Recall", "Precision", "NA", "F1 Score", "T TiTv" , "Q TiTv" , "T Het Hom" , "Q Het Hom", "Name")
melted <- melt(merged, id.vars=c("Name", "Filter", "Type"))

metrics <- subset(melted, variable%in%c("Recall", "Precision", "F1 Score"))
p1 <- ggplot(metrics, aes(x=Name, y=value, color=Filter)) +
  geom_point(stat="identity", position = position_jitter(w = 0.06, h = 0)) +
  geom_text(aes(label=ifelse(Filter=="PASS", round(value, 3), "")), color="black", size=2.5, hjust=-0.4, vjust=0.5) +
  geom_text(aes(label=ifelse(Filter!="PASS", round(value, 3), "")), color="darkgrey", size=2.5, hjust=1.6, vjust=0.5) +
  facet_grid( variable ~ Type, scales="free_y" ) +
  ylab("Metrics") +
  theme(axis.text.x=element_text(angle=30, hjust = 1))

counts <- subset(melted, variable%in%c("True Positives", "False Negatives", "False Positives"))
p2 <- ggplot(counts, aes(x=Name, y=value, color=Filter)) +
  geom_point(stat="identity", position = position_jitter(w = 0.06, h = 0))  +
  facet_grid( variable ~ Type, scales="free_y" ) +
  ylab("Counts") +
  geom_text(aes(label=ifelse(Filter=="PASS", round(value, round_digits), "")), color="black", size=2.5, hjust=-0.4, vjust=0.5) +
  geom_text(aes(label=ifelse(Filter!="PASS", round(value, round_digits), "")), color="darkgrey", size=2.5, hjust=1.6, vjust=0.5) +
  theme(axis.text.x=element_text(angle=30, hjust = 1))

ggsave(plot=multiplot(p1, p2, cols=2), filename = 'metrics.png', width=4, height=3)