#NOTE: the Java wrapper for this script first sources CNVPlottingLibrary.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

library(optparse)
option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--tangent_file", "-tangent_file"), dest="tangent_file", action="store"),
    make_option(c("--pre_tangent_file", "-pre_tangent_file"), dest="pre_tangent_file", action="store"),
    make_option(c("--segments_file", "-segments_file"), dest="segments_file", action="store"),
    make_option(c("--contig_names", "-contig_names"), dest="contig_names", action="store"),         #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--contig_lengths", "-contig_lengths"), dest="contig_lengths", action="store"),   #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--output_prefix", "-output_prefix"), dest="output_prefix", action="store"),
    make_option(c("--is_log2_input", "-is_log2_input"), dest="is_log2_input", action="store"))

opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["sample_name"]]
tangent_file = opt[["tangent_file"]]
pre_tangent_file = opt[["pre_tangent_file"]]
segments_file = opt[["segments_file"]]
contig_names_string = opt[["contig_names"]]
contig_lengths_string = opt[["contig_lengths"]]
output_dir = opt[["output_dir"]]
output_prefix = opt[["output_prefix"]]
is_log2_input = as.logical(opt[["is_log2_input"]])

#check that input files exist; if not, quit with error code that GATK will pick up
if (!all(file.exists(c(tangent_file, pre_tangent_file, segments_file)))) {
    quit(save="no", status=1, runLast=FALSE)
}

contig_names = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]])
contig_lengths = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))

CalculateQc = function(dat) {
    return(median(abs(diff(dat))))
}

#plotting is extracted to a function for debugging purposes
create_copy_ratio_plots_file = function(sample_name, tangent_file, pre_tangent_file, segments_file, contig_names, contig_lengths, output_dir, output_prefix, is_log2_input) {
    #read in file and extract needed data
    segments = read.table(segments_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    tn = read.table(tangent_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    preTn = read.table(pre_tangent_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)

    #convert the sample name field to "VALUE" for uniformity
    names(tn)[names(tn) == sample_name] = "VALUE"
    names(preTn)[names(preTn) == sample_name] = "VALUE"

    #transform tn/preTN if necessary
    if (is_log2_input) {
        tn$VALUE = 2^tn$VALUE
        preTn$VALUE = 2^preTn$VALUE
    }

    #write the QC files
    preQc = CalculateQc(preTn$VALUE)
    postQc = CalculateQc(tn$VALUE)
    write.table(round(preQc, 3), file.path(output_dir, paste(output_prefix, "_preQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(postQc, 3), file.path(output_dir, paste(output_prefix, "_postQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(preQc - postQc, 3), file.path(output_dir, paste(output_prefix, "_dQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round((preQc - postQc) / preQc, 3), file.path(output_dir, paste(output_prefix, "_scaled_dQc.txt", sep="")), col.names=FALSE, row.names=FALSE)

    #plot tangent-normalized coverage with segments
    plot_segments_file = file.path(output_dir, paste(output_prefix, "_FullGenome.png", sep=""))
    png(plot_segments_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    SetUpPlot("Tangent-normalized copy ratio", 0, 4, "", contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatioWithSegments(tn, segments, contig_names, contig_starts, FALSE)
    dev.off()

    #plot pre- and post-tangent normalization coverage on top of each other
    pre_color_blue="#3B5DFF"
    post_color_green="#4FC601"
    #plot over full range
    plot_pre_post_full_file = file.path(output_dir, paste(output_prefix, "_Before_After.png", sep=""))
    png(plot_pre_post_full_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2,1), cex=0.75, las=1)
    SetUpPlot("Copy ratio", 0, max(preTn$VALUE), paste("Pre Tangent Normalization, QC = ", round(preQc, 3), sep=""), contig_names, contig_starts, contig_ends, FALSE)
    PlotCopyRatio(preTn, pre_color_blue, contig_names, contig_starts)
    SetUpPlot("Copy ratio", 0, max(tn$VALUE), paste("Post Tangent Normalization, QC = ", round(postQc, 3), sep=""), contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatio(tn, post_color_green, contig_names, contig_starts)
    dev.off()
    #plot up to CR = 4
    plot_pre_post_CR_lim_file = file.path(output_dir, paste(output_prefix, "_Before_After_CR_Lim_4.png", sep=""))
    png(plot_pre_post_CR_lim_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2,1), cex=0.75, las=1)
    SetUpPlot("Copy ratio", 0, 4, paste("Pre Tangent Normalization, QC = ", round(preQc, 3), sep=""), contig_names, contig_starts, contig_ends, FALSE)
    PlotCopyRatio(preTn, pre_color_blue, contig_names, contig_starts)
    SetUpPlot("Copy ratio", 0, 4, paste("Post Tangent Normalization, QC = ", round(postQc, 3), sep=""), contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatio(tn, post_color_green, contig_names, contig_starts)
    dev.off()

    #check for created files and quit with error code if not found
    if (!all(file.exists(c(plot_segments_file, plot_pre_post_full_file, plot_pre_post_CR_lim_file)))) {
        quit(save="no", status=1, runLast=FALSE)
    }
}

create_copy_ratio_plots_file(sample_name, tangent_file, pre_tangent_file, segments_file, contig_names, contig_lengths, output_dir, output_prefix, is_log2_input)