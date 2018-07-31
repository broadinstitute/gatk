#NOTE: the Java wrapper for this script first sources CNVPlottingLibrary.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

library(optparse)
library(data.table)

option_list = list(
make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
make_option(c("--denoised_copy_ratios_file", "-denoised_copy_ratios_file"), dest="denoised_copy_ratios_file", action="store"),
make_option(c("--allelic_counts_file", "-allelic_counts_file"), dest="allelic_counts_file", action="store"),
make_option(c("--modeled_segments_file", "-modeled_segments_file"), dest="modeled_segments_file", action="store"),
make_option(c("--contig_names", "-contig_names"), dest="contig_names", action="store"),         #string with elements separated by "CONTIG_DELIMITER"
make_option(c("--contig_lengths", "-contig_lengths"), dest="contig_lengths", action="store"),   #string with elements separated by "CONTIG_DELIMITER"
make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
make_option(c("--output_prefix", "-output_prefix"), dest="output_prefix", action="store"))
opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["sample_name"]]
denoised_copy_ratios_file = opt[["denoised_copy_ratios_file"]]
allelic_counts_file = opt[["allelic_counts_file"]]
modeled_segments_file = opt[["modeled_segments_file"]]
contig_names_string = opt[["contig_names"]]
contig_lengths_string = opt[["contig_lengths"]]
output_dir = opt[["output_dir"]]
output_prefix = opt[["output_prefix"]]

#check that required input files exist; if not, quit with error code that GATK will pick up
if (!file.exists(modeled_segments_file)) {
    quit(save="no", status=1, runLast=FALSE)
}

contig_names = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]])
contig_lengths = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))

#plotting is extracted to a function for debugging purposes
WriteModeledSegmentsPlot = function(sample_name, allelic_counts_file, denoised_copy_ratios_file, modeled_segments_file, contig_names, contig_lengths, output_dir, output_prefix) {
    modeled_segments_df = ReadTSV(modeled_segments_file)

    plot_file = file.path(output_dir, paste(output_prefix, ".modeled.png", sep=""))
    num_plots = ifelse(all(file.exists(c(denoised_copy_ratios_file, allelic_counts_file))), 2, 1)
    png(plot_file, 12, 3.5 * num_plots, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(num_plots, 1), cex=0.75, las=1)

    if (file.exists(denoised_copy_ratios_file) && denoised_copy_ratios_file!="null") {
        denoised_copy_ratios_df = ReadTSV(denoised_copy_ratios_file)

        #transform to linear copy ratio
        denoised_copy_ratios_df[["COPY_RATIO"]] = 2^denoised_copy_ratios_df[["LOG2_COPY_RATIO"]]

        #determine copy-ratio midpoints
        denoised_copy_ratios_df[["MIDDLE"]] = round((denoised_copy_ratios_df[["START"]] + denoised_copy_ratios_df[["END"]]) / 2)

        SetUpPlot(sample_name, "denoised copy ratio", 0, 4, "contig", contig_names, contig_starts, contig_ends, TRUE)
        PlotCopyRatiosWithModeledSegments(denoised_copy_ratios_df, modeled_segments_df, contig_names, contig_starts)
    }

    if (file.exists(allelic_counts_file) && allelic_counts_file!="null") {
        allelic_counts_df = ReadTSV(allelic_counts_file)

        SetUpPlot(sample_name, "alternate-allele fraction", 0, 1.0, "contig", contig_names, contig_starts, contig_ends, TRUE)
        PlotAlternateAlleleFractionsWithModeledSegments(allelic_counts_df, modeled_segments_df, contig_names, contig_starts)
    }

    dev.off()

    #check for created file and quit with error code if not found
    if (!file.exists(plot_file)) {
        quit(save="no", status=1, runLast=FALSE)
    }
}

WriteModeledSegmentsPlot(sample_name, allelic_counts_file, denoised_copy_ratios_file, modeled_segments_file, contig_names, contig_lengths, output_dir, output_prefix)