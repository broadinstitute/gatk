#NOTE: the Java wrapper for this script first sources CNVPlottingLibrary.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

library(optparse)
option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--snp_counts_file", "-snp_counts_file"), dest="snp_counts_file", action="store"),
    make_option(c("--tangent_file", "-tangent_file"), dest="tangent_file", action="store"),
    make_option(c("--segments_file", "-segments_file"), dest="segments_file", action="store"),
    make_option(c("--contig_names", "-contig_names"), dest="contig_names", action="store"),         #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--contig_lengths", "-contig_lengths"), dest="contig_lengths", action="store"),   #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--output_prefix", "-output_prefix"), dest="output_prefix", action="store"))
opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["sample_name"]]
snp_counts_file = opt[["snp_counts_file"]]
tangent_file = opt[["tangent_file"]]
segments_file = opt[["segments_file"]]
contig_names_string = opt[["contig_names"]]
contig_lengths_string = opt[["contig_lengths"]]
output_dir = opt[["output_dir"]]
output_prefix = opt[["output_prefix"]]

#check that input files exist; if not, quit with error code that GATK will pick up
if (!all(file.exists(c(snp_counts_file, tangent_file, segments_file)))) {
    quit(save="no", status=1, runLast=FALSE)
}

contig_names = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]])
contig_lengths = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))

#plotting is extracted to a function for debugging purposes
create_acnv_plots_file = function(sample_name, snp_counts_file, tangent_file, segments_file, contig_names, contig_lengths, output_dir, output_prefix) {
    #set up coverage, snps, and segments data frames
    snp_counts = read.table(snp_counts_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    segments = read.table(segments_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    coverage = read.table(tangent_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)

    #convert the sample name field to "VALUE" for uniformity
    names(coverage)[names(coverage) == sample_name] = "VALUE"

    #ACNV coverage is always in log2 space
    coverage$VALUE = 2^coverage$VALUE

    #plot ACNV CR and MAF data and segment posteriors
    plot_file = file.path(output_dir, paste(output_prefix, "_ACNV.png", sep=""))
    png(plot_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2,1), cex=0.75, las=1)
    SetUpPlot("Tangent-normalized copy ratio", 0, 4, "Contig", contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatioWithSegments(coverage, segments, contig_names, contig_starts, TRUE)
    SetUpPlot("Minor-allele fraction", 0, 0.5, "Contig", contig_names, contig_starts, contig_ends, TRUE)
    PlotAlleleFractionWithSegments(snp_counts, segments, contig_names, contig_starts)
    dev.off()

    #check for created file and quit with error code if not found
    if (!file.exists(plot_file)) {
        quit(save="no", status=1, runLast=FALSE)
    }
}

create_acnv_plots_file(sample_name, snp_counts_file, tangent_file, segments_file, contig_names, contig_lengths, output_dir, output_prefix)

