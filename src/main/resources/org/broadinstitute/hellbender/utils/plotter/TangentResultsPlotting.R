#NOTE: the Java wrapper for this script first sources plotting.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

library(optparse)
option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--targets_file", "-targets_file"), dest="targets_file", action="store"),
    make_option(c("--pre_tn_file", "-pre_tn_file"), dest="pre_tn_file", action="store"),
    make_option(c("--segments_file", "-segments_file"), dest="segments_file", action="store"),
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--output_prefix", "-output_prefix"), dest="output_prefix", action="store"),
    make_option(c("--log2_input", "-log"), dest="log2_input", action="store"),
    make_option(c("--sex_chrs", "-sexchrs"), dest="sex_chrs", action="store"))

opt = parse_args(OptionParser(option_list=option_list))
save(opt, file="debug.RData")

sample_name=opt[["sample_name"]]
tn_file=opt[["targets_file"]]
pre_tn_file=opt[["pre_tn_file"]]
segments_file=opt[["segments_file"]]
output_dir=opt[["output_dir"]]
output_prefix=opt[["output_prefix"]]
log_input=as.logical(opt[["log2_input"]])
sex_chrs=as.logical(opt[["sex_chrs"]])
num_chromosomes = ifelse(sex_chrs, 24, 22)

#check input files exist.  If not, quit with error code that GATK will pick up
if (!all(file.exists(c(tn_file, pre_tn_file, segments_file)))) {
    quit(save = "no", status = 1, runLast = FALSE)
}

QC = function(dat) {
    return(median(abs(diff(dat))))
}

# Use a function for debugging purposes
create_tangent_plots_file = function(sample_name, tn_file, pre_tn_file, segments_file, output_dir, output_prefix, log_input, sex_chrs) {
	# Read in file and extract needed data
	segments = read.table(segments_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	tn = read.table(tn_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	preTn = read.table(pre_tn_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	names(tn)[names(tn) == sample_name] = "VALUE"
	names(preTn)[names(preTn) == sample_name] = "VALUE"

	if (log_input) {
	    tn$VALUE = 2^tn$VALUE
	    preTn$VALUE = 2^preTn$VALUE
	}

    preQc = QC(preTn$VALUE)
    postQc = QC(tn$VALUE)

    write.table(round(preQc, 3), file.path( output_dir, paste(output_prefix, "_preQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(postQc, 3), file.path( output_dir, paste(output_prefix, "_postQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(preQc - postQc, 3), file.path( output_dir, paste(output_prefix, "_dQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round((preQc - postQc)/preQc, 3), file.path( output_dir, paste(output_prefix, "_scaled_dQc.txt", sep="")), col.names=FALSE, row.names=FALSE)

    #plot tangent-normalized coverage with segments
    plot.fn = file.path( output_dir, paste(output_prefix, "_FullGenome.png", sep=""))
    png(plot.fn, 12, 7, units="in", type="cairo", res=300, bg="white")
    SetUpPlot("Total copy ratio", 0, 4, "", TRUE, num_chromosomes)
    PlotCopyRatioWithSegments(tn, segments, FALSE)
    dev.off()

    pre_color_blue="#3B5DFF"
    post_color_green="#4FC601"

    #plot pre- and post-tangent normalization coverage on top of each other
    plot.fns = c(file.path(output_dir, paste(output_prefix, "_Before_After_CR_Lim_4.png", sep="")),
            file.path(output_dir, paste(output_prefix, "_Before_After.png", sep="")))

    for (i in 1:2) {
        png(plot.fns[i], 12, 7, units="in", type="cairo", res=300, bg="white")
        par(mfrow=c(2,1), cex=0.75, las=1)
        SetUpPlot("Copy-Ratio", 0, ifelse(i == 1, 4, max(preTn$VALUE)), paste("Pre Tangent Normalization, QC = ", round(preQc, 3), sep=""), FALSE, num_chromosomes)
        PlotCopyRatio(preTn, pre_color_blue)
        SetUpPlot("Copy-Ratio", 0, ifelse(i == 1, 4, max(tn$VALUE)), paste("Post Tangent Normalization, QC = ", round(postQc, 3), sep=""), TRUE, num_chromosomes)
        PlotCopyRatio(tn, post_color_green)
        dev.off()
    }

    #check for created files and quit with error code if not found
    if (!all(file.exists(c(plot.fns, plot.fn)))) {
       quit(save = "no", status = 1, runLast = FALSE)
    }
}

create_tangent_plots_file(sample_name, tn_file, pre_tn_file, segments_file, output_dir, output_prefix, log_input, sex_chrs)