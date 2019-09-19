library(dplyr)
library(ggplot2)

inverseLogitFn = function(score) {
    v = 10.0^(-score/10.0)
    v / (1.0 + v)
}

args <- commandArgs(TRUE)
verbose = TRUE

variant_eval_filename <- args[1]
output_basename <- tools::file_path_sans_ext(args[2])
sample <- args[3]


make_sample_df = function(filename, sample) {
    df = read.table(filename, sep="", stringsAsFactors = F, header = T)
    df = df %>% rowwise() %>% dplyr::filter(Filter == "called" ) # only called sites
    df = df %>% rowwise() %>% 
	    transmute(AF_bin = inverseLogitFn(AlleleFrequency),
                      EvaluationType = EvalFeatureInput,
                      avgVarAF= avgVarAF)

    merge(subset(df, EvaluationType == "eval", select = -EvaluationType),
          subset(df, EvaluationType == "thousand_genomes", select = -EvaluationType),
          by = c("AF_bin"), suffixes = c(".array", ".thousand_g")) %>%
      mutate(sample = sample)
}

run_single_sample = function(filename, output_basename, sample) {
    single_samp = make_sample_df(filename, sample)

    png(paste(output_basename,  ".af_differences", ".png", sep = ""), width = 6, height = 6, units = 'in', res = 300)
    print(ggplot(single_samp, aes(x = avgVarAF.thousand_g, y = avgVarAF.array)) + geom_point() +
        geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.7) +
        labs(x = "AF from Thousand Genomes", y = "AF from Array"))
    dev.off()

    png(paste(output_basename,  ".af", ".png", sep = ""), width = 6, height = 6, units = 'in', res = 300)
    print(ggplot(single_samp) + geom_point(aes(x = AF_bin, y = avgVarAF.array, color = "array")) +
        geom_point(aes(x = AF_bin, y = avgVarAF.thousand_g, color = "thousand genomes"))  +
        geom_abline(slope = 1, intercept = 0, color = "black", alpha = 0.7) +
        labs(x = "AF Bin", y = "AF in VCFs") + theme(legend.position = "bottom"))
    dev.off()
}


run_single_sample(variant_eval_filename, output_basename, sample)
