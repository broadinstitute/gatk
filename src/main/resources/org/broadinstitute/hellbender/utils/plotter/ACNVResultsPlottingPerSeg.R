#NOTE: the Java wrapper for this script first sources CNV_plotting_library.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

library(optparse)
option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--snp_counts_file", "-snp_counts_file"), dest="snp_counts_file", action="store"),
    make_option(c("--coverage_file", "-coverage_file"), dest="coverage_file", action="store"),
    make_option(c("--segments_file", "-segments_file"), dest="segments_file", action="store"),
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--output_prefix", "-output_prefix"), dest="output_prefix", action="store"),
    make_option(c("--sex_chrs", "-sexchrs"), dest="sex_chrs", action="store"))

opt = parse_args(OptionParser(option_list=option_list))
save(opt, file="debug.RData")

sample_name=opt[["sample_name"]]
snp_counts_file=opt[["snp_counts_file"]]
coverage_file=opt[["coverage_file"]]
segments_file=opt[["segments_file"]]
output_dir=opt[["output_dir"]]
output_prefix=opt[["output_prefix"]]
sex_chrs=as.logical(opt[["sex_chrs"]])
num_chromosomes = ifelse(sex_chrs, 24, 22)

#check that input files exist. If not, quit with an error code that GATK will pick up
if (!all(file.exists(c(snp_counts_file, coverage_file, segments_file)))) {
    quit(save = "no", status = 1, runLast = FALSE)
}

create_acnv_plots_file = function(sample_name, snp_counts_file, coverage_file, segments_file, output_dir, output_prefix, num_chromosomes) {
    #set up coverage, snps, and segments data frames
    snp_counts = read.table(snp_counts_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    segments = read.table(segments_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    snp_counts[, "CONTIG"] = convert_XY_to_23_24(snp_counts[, "CONTIG"])
    coverage = read.table(coverage_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)

    #convert the sample name field to "VALUE" for uniformity
    headers = names(coverage)
    headers[!headers %in% c("contig", "start", "stop", "name")] = "VALUE"
    names(coverage) = headers
    coverage$VALUE = 2^coverage$VALUE #ACNV is in log space
    segments[, "Chromosome"] = convert_XY_to_23_24(segments[, "Chromosome"])
    num_segs = length(which(segments[,"Chromosome"] <= num_chromosomes))

    #subset to autosomes if sex_chrs not specified
    is_single_seg = function(chr) sum(segments$Chromosome == chr) == 1
    single_seg_chrs = Filter(is_single_seg, c(1:num_chromosomes))

    #if the last segment is on its own segment, we need to plot it additionally
    if(length(single_seg_chrs)>0 && tail(single_seg_chrs, 1) == segments[nrow(segments), "Chromosome"]){ num_segs=num_segs+1 }
    for( i in 1:(num_segs-1) ) {
        seg1=segments[i,]
        chr = seg1[, "Chromosome"]
        only_1_seg = chr %in% single_seg_chrs
        if(! only_1_seg) {
            seg2=segments[i+1,]
            if(chr != seg2[, "Chromosome"]) {next()}
            segs = rbind(seg1, seg2)
        } else {
            segs=seg1
        }

        #make the plots
        plot_file_name = file.path(output_dir, paste(output_prefix, "_Chr_", chr, "_", segs[1, "Start"], "_",
            segs[nrow(segs), "End"], ".png", sep=""))
        png(plot_file_name, 9, 5, units="in", type="cairo", res=300, bg="white")
        par(mfrow=c(2,1), cex=0.75, las=1)

        #We need to look at the coverage to scale the y-axis
        segment_start = segs[1, "Start"]
        segment_end = segs[nrow(segs), "End"]
        target_indices = which(coverage[, "contig"] == chr & coverage[, "stop"] > segment_start &
            coverage[, "stop"] < segment_end)
        coverages = coverage[target_indices, "VALUE"]

        SetUpPlotPerSeg("Tangent-Normalized Coverage", min(coverages), max(coverages), paste("Chromosome ", chr, "   ||    num-targets/snps: ", segs[1, "Num_Targets"], "/", segs[1, "Num_SNPs"], ", ", ifelse(only_1_seg, 0, segs[2, "Num_Targets"]), "/", ifelse(only_1_seg, 0, segs[2, "Num_SNPs"]), sep=""), TRUE, FALSE, segs)
        PlotCopyRatioWithSegments(coverage, segs, TRUE, FALSE, 0.3)
        SetUpPlotPerSeg("Minor Allele Fraction", 0, 0.5, "Position (Mbp)", TRUE, TRUE, segs)
        PlotAlleleFraction(snp_counts, segs, FALSE, 0.3)
        dev.off()
    }

    #check for created file and quit with error code if not found
    if (!file.exists(plot_file_name)) {
        quit(save = "no", status = 1, runLast = FALSE)
    }
}

create_acnv_plots_file(sample_name, snp_counts_file, coverage_file, segments_file, output_dir, output_prefix, num_chromosomes)

