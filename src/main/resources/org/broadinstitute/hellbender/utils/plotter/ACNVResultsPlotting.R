#NOTE: the Java wrapper for this script first sources CNV_plotting_library.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q()}))    # Useful for debugging

library(optparse)
option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--snp_counts_file", "-snp_counts_file"), dest="snp_counts_file", action="store"),
    make_option(c("--coverage_file", "-coverage_file"), dest="coverage_file", action="store"),
    make_option(c("--segments_file", "-segments_file"), dest="segments_file", action="store"),
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--sex_chrs", "-sexchrs"), dest="sex_chrs", action="store"))

opt = parse_args(OptionParser(option_list=option_list))
save(opt, file="debug.RData")

sample_name=opt[["sample_name"]]
snp_counts_file=opt[["snp_counts_file"]]
coverage_file=opt[["coverage_file"]]
segments_file=opt[["segments_file"]]
output_file=opt[["output_dir"]]
sex_chrs=as.logical(opt[["sex_chrs"]])
num_chromosomes = ifelse(sex_chrs, 24, 22)

create_acnv_plots_file = function(sample_name, snp_counts_file, coverage_file, segments_file, output_dir, num_chromosomes) {
	#set up coverage, snps, and segments data frames
	snp_counts = read.table(snp_counts_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	segments = read.table(segments_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    snp_counts[,"CONTIG"] = convert_XY_to_23_24(snp_counts[,"CONTIG"])
    coverage = read.table(coverage_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    names(coverage)[names(coverage) == sample_name] = "VALUE"
    coverage$VALUE=2^coverage$VALUE #ACNV is in log space
    segments[,"Chromosome"] = convert_XY_to_23_24(segments[,"Chromosome"])

    #make the plots
    plot_file_name = file.path(output_dir, paste(sample_name, "_ACNV.png", sep=""))
    png(plot_file_name, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2,1), cex=0.75, las=1)
    SetUpPlot("Minor Allele Fraction", 0, 0.5, "Chromosome", TRUE, num_chromosomes)
    PlotAlleleFraction(snp_counts, segments)
    SetUpPlot("Tangent-Normalized Coverage", 0, 4, "Chromosome", TRUE, num_chromosomes)
    PlotCopyRatioWithSegments(coverage, segments, TRUE)
    dev.off()
}

create_acnv_plots_file(sample_name, snp_counts_file, coverage_file, segments_file, output_file, num_chromosomes)