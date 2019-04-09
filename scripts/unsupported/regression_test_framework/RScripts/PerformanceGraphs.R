library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(tidyr)
library(gridExtra)

loadTable <- function(input){
  newMD <- read.csv(input, stringsAsFactors = FALSE)
  return(newMD)
}

reformatData <- function(df) {
  mod <- df %>% mutate(persisient = str_replace(persisient, "p", "Preemptible")) %>% mutate(persisient = str_replace(persisient, "x", "non-Preemptible"))
  mod <- mod %>% mutate(sku_description = str_replace(sku_description, " attached to Preemptible VMs", ""))  %>% mutate(sku_description = str_replace(sku_description, "Preemptible ", ""))
  mod <- mod %>% mutate(sku_description = str_replace(sku_description, " from Americas to Americas", ""))
  mod <- mod %>% mutate(sku_description = str_replace(sku_description, " running in Americas", ""))
  mod$Preemptible <- mod$persisient
  mod$persisient <- NULL
  return(mod)
}

################################## Inputting Data ########################################


trialTotal <- "/Users/emeryj/hellbender/Scripts/HaplotypeCallerSpark/Analysis/CheckpointTotalReductionPlots.csv"
trialstreambam <- "/Users/emeryj/hellbender/Scripts/HaplotypeCallerSpark/Analysis/CheckpointRuntimeTrial3.csv"
trialrepeated <- "/Users/emeryj/hellbender/Scripts/HaplotypeCallerSpark/Analysis/RuntimeComparisonRepettions.csv"

tbl1 <- loadTable(trialTotal)
tbl3 <- loadTable(trialstreambam)

tblr <- loadTable(trialrepeated)

filtered1 <- tbl1 %>% filter(grepl("ome",Type, fixed=TRUE))
filtered1 <- filtered1 %>% mutate(Ratio = -Ratio)
filtered3 <- tbl3 %>% filter(grepl("ome",Type, fixed=TRUE))
filtered3 <- filtered3 %>% mutate(Ratio = -Ratio)

tblr <- tblr %>% filter(grepl("a",Type, fixed=TRUE))
tblr <- tblr %>% mutate(Runtime = Runtime / 60)
gatk <- tblr %>% filter(grepl("GATK",Type, fixed=TRUE))
after <- tblr %>% filter(grepl("After",Type, fixed=TRUE))
mean(after, Runtime)

p <- ggplot(data=tblr, aes(x=Runtime, fill=Type, alpha = 0.5))
p <- p + labs(x = "Total HaplotypeCaller Chr1 runtime (minutes)", y = "Number of Trials", title = "20 Repetitions of G96830.NA12878 Runtime Histogram")
p <- p + geom_histogram(binwidth = 5, position = "identity")
p <- p + geom_vline(xintercept =mean(gatk$Runtime)) + geom_vline(xintercept =mean(after$Runtime)) 
print(p)

p <- ggplot(data=tblr, aes(x=Type, y=Runtime))
p <- p + labs(x = "Gatk Jar", y = "Total HaplotypeCaller Chr1 runtime (minutes)", title = "20 Repetitions of G96830.NA12878 Runtime Histogram")
p <- p + geom_boxplot()
p <- p + coord_flip()
p <- p + theme_bw()
p <- p + scale_y_continuous(limits=c(0,230))
p <- p + geom_vline(xintercept =mean(gatk$Runtime)) + geom_vline(xintercept =mean(after$Runtime)) 
print(p)

sum(exomes1$Baseline)
median(exomes1$Ratio)
mean(filtered1$Ratio)
p <- ggplot(data=filtered1, aes(x=Title, y=Ratio))
p <- p + geom_bar(stat="identity", aes(fill=Type)) + scale_fill_grey()
p <- p + labs(x = "Trial Bam", y = "Percent change in runtime")
p <- p + theme_bw()
p <- p + theme(axis.text.x=element_text(angle=-20,hjust=1))
p <- p + scale_x_discrete(position = "top") 
p <- p + geom_hline(yintercept =mean(filtered1$Ratio))
print(p)


p <- ggplot(data=filtered3, aes(x=Title, y=Ratio))
p <- p + geom_bar(stat="identity", aes(fill=Type)) + scale_fill_grey()
p <- p + labs(x = "Trial Bam", y = "Percent change in runtime")
p <- p + theme_bw()
p <- p + theme(axis.text.x=element_text(angle=-20,hjust=1))
p <- p + scale_x_discrete(position = "top") 
p <- p + geom_hline(yintercept =mean(filtered1$Ratio))
print(p)


