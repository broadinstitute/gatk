#NOTE: the Java wrapper for this script first sources CNVPlottingLibrary.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

library(optparse)
library(data.table)

option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--standardized_copy_ratios_file", "-standardized_copy_ratios_file"), dest="standardized_copy_ratios_file", action="store"),
    make_option(c("--denoised_copy_ratios_file", "-denoised_copy_ratios_file"), dest="denoised_copy_ratios_file", action="store"),
    make_option(c("--contig_names", "-contig_names"), dest="contig_names", action="store"),         #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--contig_lengths", "-contig_lengths"), dest="contig_lengths", action="store"),   #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--output_prefix", "-output_prefix"), dest="output_prefix", action="store"))

opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["sample_name"]]
standardized_copy_ratios_file = opt[["standardized_copy_ratios_file"]]
denoised_copy_ratios_file = opt[["denoised_copy_ratios_file"]]
contig_names_string = opt[["contig_names"]]
contig_lengths_string = opt[["contig_lengths"]]
output_dir = opt[["output_dir"]]
output_prefix = opt[["output_prefix"]]

#check that input files exist; if not, quit with error code that GATK will pick up
if (!all(file.exists(c(standardized_copy_ratios_file, denoised_copy_ratios_file)))) {
    quit(save="no", status=1, runLast=FALSE)
}

contig_names = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]])
contig_lengths = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))

CalculateMedianAbsoluteDeviation = function(dat) {
    return(median(abs(diff(dat))))
}

#plotting is extracted to a function for debugging purposes
WriteDenoisingPlots = function(sample_name, standardized_copy_ratios_file, denoised_copy_ratios_file, contig_names, output_dir, output_prefix) {
    standardized_copy_ratios_df = ReadTSV(standardized_copy_ratios_file)
    denoised_copy_ratios_df = ReadTSV(denoised_copy_ratios_file)

    #transform to linear copy ratio
    standardized_copy_ratios_df[["COPY_RATIO"]] = 2^standardized_copy_ratios_df[["LOG2_COPY_RATIO"]]
    denoised_copy_ratios_df[["COPY_RATIO"]] = 2^denoised_copy_ratios_df[["LOG2_COPY_RATIO"]]

    #determine copy-ratio midpoints
    standardized_copy_ratios_df[["MIDDLE"]] = round((standardized_copy_ratios_df[["START"]] + standardized_copy_ratios_df[["END"]]) / 2)
    denoised_copy_ratios_df[["MIDDLE"]] = round((denoised_copy_ratios_df[["START"]] + denoised_copy_ratios_df[["END"]]) / 2)

    #write the MAD files
    standardizedMAD = CalculateMedianAbsoluteDeviation(standardized_copy_ratios_df[["COPY_RATIO"]])
    denoisedMAD = CalculateMedianAbsoluteDeviation(denoised_copy_ratios_df[["COPY_RATIO"]])
    write.table(round(standardizedMAD, 3), file.path(output_dir, paste(output_prefix, ".standardizedMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(denoisedMAD, 3), file.path(output_dir, paste(output_prefix, ".denoisedMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(standardizedMAD - denoisedMAD, 3), file.path(output_dir, paste(output_prefix, ".deltaMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round((standardizedMAD - denoisedMAD) / standardizedMAD, 3), file.path(output_dir, paste(output_prefix, ".scaledDeltaMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)

    #plot standardized and denoised copy ratio on top of each other
    pre_color_blue = "#3B5DFF"
    post_color_green = "#4FC601"

    #plot over full range
    denoising_plot_file = file.path(output_dir, paste(output_prefix, ".denoised.png", sep=""))
    png(denoising_plot_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2, 1), cex=0.75, las=1)
    SetUpPlot(sample_name, "standardized copy ratio", 0, max(standardized_copy_ratios_df[["COPY_RATIO"]]), paste("median absolute deviation = ", round(standardizedMAD, 3), sep=""), contig_names, contig_starts, contig_ends, FALSE)
    PlotCopyRatios(standardized_copy_ratios_df, pre_color_blue, contig_names, contig_starts)
    SetUpPlot(sample_name, "denoised copy ratio", 0, max(denoised_copy_ratios_df[["COPY_RATIO"]]), paste("median absolute deviation = ", round(denoisedMAD, 3), sep=""), contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatios(denoised_copy_ratios_df, post_color_green, contig_names, contig_starts)
    dev.off()

    #plot up to CR = 4
    denoising_limit_plot_file = file.path(output_dir, paste(output_prefix, ".denoisedLimit4.png", sep=""))
    png(denoising_limit_plot_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2, 1), cex=0.75, las=1)
    SetUpPlot(sample_name, "standardized copy ratio", 0, 4, paste("median absolute deviation = ", round(standardizedMAD, 3), sep=""), contig_names, contig_starts, contig_ends, FALSE)
    PlotCopyRatios(standardized_copy_ratios_df, pre_color_blue, contig_names, contig_starts)
    SetUpPlot(sample_name, "denoised copy ratio", 0, 4, paste("median absolute deviation = ", round(denoisedMAD, 3), sep=""), contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatios(denoised_copy_ratios_df, post_color_green, contig_names, contig_starts)
    dev.off()

    #check for created files and quit with error code if not found
    if (!all(file.exists(c(denoising_plot_file, denoising_limit_plot_file)))) {
        quit(save="no", status=1, runLast=FALSE)
    }
}

WriteDenoisingPlots(sample_name, standardized_copy_ratios_file, denoised_copy_ratios_file, contig_names, output_dir, output_prefix)