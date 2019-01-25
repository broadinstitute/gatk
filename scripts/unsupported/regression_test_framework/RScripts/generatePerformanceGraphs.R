# Title     : TODO
# Objective : TODO
# Created by: emeryj
# Created on: 2/5/19

install.packages("ggplot2", repos="http://cran.us.r-project.org")
install.packages("reshape2", repos="http://cran.us.r-project.org")
install.packages("dplyr", repos="http://cran.us.r-project.org")

library(ggplot2)
library(reshape2)
library(dplyr)

loadTable <- function(input){
    newMD <- read.csv(input, stringsAsFactors = FALSE)
    return(newMD)
}

bind_save_function <- function(outputdir, prefix){
    save_function <- function(name,height=10, width=10, plot=last_plot()) {
        name_pieces <- c(outputdir, "/",prefix,"",name, ".pdf")
        filename <- paste(name_pieces, collapse='')

        print(paste("Saving",filename))
        ggsave(file=filename, height=height, width=width, units="in", plot=plot, limitsize=FALSE)
    }
}

args <- commandArgs(trailingOnly = TRUE)
print(args)

input_csv = args[1]
input_bam = args[2]

tblr <- loadTable(input_csv)

tblr <- tblr %>% filter(grepl("a",Type, fixed=TRUE))
tblr <- tblr %>% mutate(Runtime = Runtime / 60)
gatk <- tblr %>% filter(grepl("GATK",Type, fixed=TRUE)) ## todo this will need updated values
after <- tblr %>% filter(grepl("After",Type, fixed=TRUE))

save_with_name <- bind_save_function(".", "plot")

p <- ggplot(data=tblr, aes(x=Runtime, fill=Type, alpha = 0.5))
p <- p + labs(x = "Total HaplotypeCaller Chr1 runtime (minutes)", y = "Number of Trials", title = paste("20 Repetitions of ",input_bam," Runtime Histogram"))
p <- p + geom_histogram(binwidth = 5, position = "identity")
p <- p + geom_vline(xintercept =mean(gatk$Runtime)) + geom_vline(xintercept =mean(after$Runtime))
print(p)

save_with_name("hisogram")

p <- ggplot(data=tblr, aes(x=Type, y=Runtime))
p <- p + labs(x = "Gatk Jar", y = "Total HaplotypeCaller Chr1 runtime (minutes)", title = paste("20 Repetitions of ",input_bam," Runtime Boxplot"))
p <- p + geom_boxplot()
p <- p + coord_flip()
p <- p + theme_bw()
p <- p + scale_y_continuous(limits=c(0,mean(after, Runtime)+10))
p <- p + geom_vline(xintercept =mean(gatk$Runtime)) + geom_vline(xintercept =mean(after$Runtime))
print(p)
save_with_name("boxplot")