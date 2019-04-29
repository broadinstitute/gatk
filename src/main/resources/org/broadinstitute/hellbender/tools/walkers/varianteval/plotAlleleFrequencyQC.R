library(plyr)
library(dplyr)
library(ggplot2)

inverseLogitFn = function(score) {
    v = 10.0^(-score/10.0)
    v / (1.0 + v)
}

probToPhred = function(prob){
    as.integer(- 10. * log10(prob))
}

args <- commandArgs(TRUE)
verbose = TRUE

filename <- args[1]
sample <- args[3]

folder <- gsub("(.*)\\/.*", "\\1", args[2]) ## just the path / directory
output_basename <- paste(folder, "/", sample, sep="")

run_single_sample = function(filename, output_basename, sample) {

    df = read.table(filename, sep = "", stringsAsFactors = F, header = T) %>%
        filter(Filter == "called" ) %>% rowwise()  %>%
        mutate(AF = inverseLogitFn(AlleleFrequency),
        num_sites = totalHetSites + totalHomRefSites + totalHomVarSites,
        EvalFeatureInput = ifelse(EvalFeatureInput == "eval", "Array VCF", "Thousand Genomes VCF"))

    test_results = ks.test(subset(df, EvalFeatureInput == "Array VCF")$avgVarAlleles,
    subset(df, EvalFeatureInput == "Thousand Genomes VCF")$avgVarAlleles)
    output_test_df = data.frame(SAMPLE = sample, TEST_TYPE = "Variant Allele Frequency",
    `LOD SCORE` = probToPhred(test_results$p.value),
    PVALUE = round(test_results$p.value, digits = 4),
    `FILTER STATUS` = "Called Sites Only")

    png(filename = paste(output_basename, ".af_differences.png", sep = ""), width = 6, height = 6, units = 'in', res = 300)
    print(ggplot(df, aes(x = AF, y = avgVarAlleles, color = EvalFeatureInput, size = totalCalledAlleles)) +
        geom_point() + theme(legend.position = "bottom") +
        labs(x = "Expected AF from Thousand Genomes", y = "Observed AF", title = "Expected vs Actual Allele Frequencies"))
    dev.off()

    png(filename = paste(output_basename,".af_cdf.png", sep = ""), width = 6, height = 6, units = 'in', res = 300)
    print(ggplot(df, aes(avgVarAlleles, color = EvalFeatureInput)) + stat_ecdf(alpha=0.7) + theme(legend.position = "bottom") +
        labs(x = "Observed AF", y = "Density", title = "CDF of Allele Frequencies: Thousand Genomes vs Array VCF"))
    dev.off()

    df_with_calculated_vals = df %>% rowwise() %>% transmute(AF = AF, EvalFeatureInput = EvalFeatureInput,
    p_2 = totalHomVarSites/num_sites,
    two_pq = totalHetSites/num_sites,
    q_2 = totalHomRefSites/num_sites)

    png(filename = paste(output_basename, ".p^2_differences.png", sep = ""), width = 6, height = 6, units = 'in', res = 300)
    print(ggplot(df_with_calculated_vals, aes(x=AF, y=p_2, color = EvalFeatureInput)) + geom_point() + theme(legend.position = "bottom") +
        labs(x = "AF bin from thousand genomes", y = "Fraction of sites that are homozygous variant(p^2)"))
    dev.off()

    png(filename = paste(output_basename, ".2pq_differences.png", sep = ""), width = 6, height = 6, units = 'in', res = 300)
    print(ggplot(df_with_calculated_vals, aes(x = AF, y = two_pq, color = EvalFeatureInput)) + geom_point() + theme(legend.position = "bottom") +
        labs(x = "AF bin from thousand genomes", y = "Fraction of sites that are heterozygous (2pq)"))

    dev.off()

    heterozygosity_p_2 = ks.test(subset(df_with_calculated_vals, EvalFeatureInput == "Array VCF")$p_2,
    subset(df_with_calculated_vals, EvalFeatureInput == "Thousand Genomes VCF")$p_2)
    heterozygosity_2pq = ks.test(subset(df_with_calculated_vals, EvalFeatureInput == "Array VCF")$two_pq,
    subset(df_with_calculated_vals, EvalFeatureInput == "Thousand Genomes VCF")$two_pq)

    output_test_df = rbind(output_test_df,
    data.frame(SAMPLE = sample, TEST_TYPE = "Excess Variant Homozygosity (p^2)", `LOD SCORE` = probToPhred(heterozygosity_p_2$p.value), PVALUE=round(heterozygosity_p_2$p.value, digits=4), `FILTER STATUS`="Called Sites Only"),
    data.frame(SAMPLE = sample, TEST_TYPE = "Excess Heterozygosity (2pq)", `LOD SCORE` = probToPhred(heterozygosity_2pq$p.value), PVALUE=round(heterozygosity_2pq$p.value, digits=4), `FILTER STATUS`="Called Sites Only"))

    write.table(output_test_df, sep = "\t", file = paste(output_basename, ".pvals.txt", sep = ""), row.names = F, quote = F)

    # make sure tests fail if significant p value
    if (any(output_test_df$PVALUE <= 0.05)) {
        significant_pvals =  paste(unlist(as.character(subset(output_test_df, PVALUE <= 0.05)$TEST_TYPE)), collapse = " and ")
        stop(paste("The p-value from the Kolmogorov-Smirnov test is less than 0.05: the array VCF has significantly different ", significant_pvals, " compared to the given Thousand Genomes VCF.", sep=""))
    }
}

run_single_sample(filename, output_basename, sample)