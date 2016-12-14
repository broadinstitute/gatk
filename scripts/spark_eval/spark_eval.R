if(!require(dplyr)) {
  install.packages("dplyr")
}
if(!require(ggplot2)) {
  install.packages("ggplot2")
}
library(dplyr)
library(ggplot2)

# Test Case #2: CountReadsSpark, no intervals, HDFS input
test_case_2 <- read.csv("/Users/tom/workspace/gatk/q4_spark_eval/test_case_2_results.csv")
test_case_2$Total.cores <- with(test_case_2, Number.of.executors * Executor.cores)
test_case_2$cm <- with(test_case_2, Total.cores / Time..mins.)
test_case_2$Executor.cores <- as.factor(test_case_2$Executor.cores)

png("/Users/tom/workspace/gatk/q4_spark_eval/test_case_2_time.png")
qplot(Total.cores, Time..mins., data = test_case_2, geom = c("point", "line"), colour = Executor.cores,
      xlab = "Cores", ylab = "Time (mins)",
      main = "Test Case #2: CountReadsSpark, no intervals, HDFS input")
dev.off()

speedup <- test_case_2 %>%
  filter(Total.cores == 4) %>%
  summarize(mean4cores = mean(Time..mins.)) %>%
  as.numeric()
test_case_2$Speedup <- with(test_case_2, 4 * speedup / Time..mins.)
png("/Users/tom/workspace/gatk/q4_spark_eval/test_case_2_speedup.png")
qplot(Total.cores, Speedup, data = test_case_2, geom = c("point", "line"), colour = Executor.cores,
      xlab = "Cores", ylab = "Speedup",
      main = "Test Case #2: CountReadsSpark, no intervals, HDFS input")
dev.off()

# Test Case #4: CountReadsSpark, with intervals, HDFS input
test_case_4 <- read.csv("/Users/tom/workspace/gatk/q4_spark_eval/test_case_4_results.csv")
test_case_4$Total.cores <- with(test_case_4, Number.of.executors * Executor.cores)
test_case_4$cm <- with(test_case_4, Total.cores / Time..mins.)
test_case_4$Executor.cores <- as.factor(test_case_4$Executor.cores)

png("/Users/tom/workspace/gatk/q4_spark_eval/test_case_4_time.png")
qplot(Total.cores, Time..mins., data = test_case_4, geom = c("point", "line"), colour = Executor.cores,
      xlab = "Cores", ylab = "Time (mins)",
      main = "Test Case #4: CountReadsSpark, with intervals, HDFS input")
dev.off()

speedup <- test_case_4 %>%
  filter(Total.cores == 4) %>%
  summarize(mean4cores = mean(Time..mins.)) %>%
  as.numeric()
test_case_4$Speedup <- with(test_case_4, 4 * speedup / Time..mins.)
png("/Users/tom/workspace/gatk/q4_spark_eval/test_case_4_speedup.png")
qplot(Total.cores, Speedup, data = test_case_4, geom = c("point", "line"), colour = Executor.cores,
      xlab = "Cores", ylab = "Speedup",
      main = "Test Case #4: CountReadsSpark, with intervals, HDFS input")
dev.off()

