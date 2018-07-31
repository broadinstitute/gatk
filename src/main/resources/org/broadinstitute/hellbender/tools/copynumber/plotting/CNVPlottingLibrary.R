ReadTSV = function(tsv_file) {
    # We need to filter out header lines beginning with '@';
    # however, the standard 'fread("grep ...")' causes issues with the default Docker container, so we use a temporary file.
    # See https://github.com/broadinstitute/gatk/issues/4140.
    temp_file = tempfile()
    system(sprintf("grep -v ^@ %s > %s", tsv_file, temp_file))
    return(suppressWarnings(fread(temp_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE)))
}

## Contig Background for data
SetUpPlot = function(sample_name, y.lab, y.min, y.max, x.lab, contig_names, contig_starts, contig_ends, do_label_contigs) {
    num_contigs = length(contig_names)
    contig_centers = (contig_starts + contig_ends) / 2
    genome_length = contig_ends[num_contigs]
    suppressWarnings(par(mar=c(3.1, 3.6, 3.6, 0), mgp=c(2, -0.2, -1.1)))
    plot(0, type="n", bty="n", xlim=c(0, genome_length), ylim=c(0, y.max), xlab="", ylab="", main=sample_name, xaxt="n")
    mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
    mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)

    if (do_label_contigs) {
        mtext(text=contig_names[1:num_contigs], side=1, line=ifelse(c(1:num_contigs) %% 2 == 1, -0.45, 0.0),
        at = contig_centers[1:num_contigs], las=1, cex=par("cex.axis") * par("cex") * 0.7)
    }

    for (i in 1:num_contigs) {
        use.col = ifelse(i %% 2 == 1, "grey90", "white")
        rect(xleft=contig_starts[i], ybottom=y.min, xright=contig_ends[i], ytop=y.max, col=use.col, border=NA)
    }
}

PlotCopyRatios = function(copy_ratios_df, color, contig_names, contig_starts) {
    genomic_coordinates = contig_starts[match(copy_ratios_df[["CONTIG"]], contig_names)] + copy_ratios_df[["MIDDLE"]]
    points(x=genomic_coordinates, y=copy_ratios_df[["COPY_RATIO"]], col=color, pch=".", cex=0.2)
}

PlotCopyRatiosWithModeledSegments = function(denoised_copy_ratios_df, modeled_segments_df, contig_names, contig_starts, point_size=0.2) {
   points_start_index = 1
   for (s in 1:nrow(modeled_segments_df)) {
       #skip segments with no points
       num_points = modeled_segments_df[s, "NUM_POINTS_COPY_RATIO"]
       if (num_points == 0) {
           next
       }
       points_end_index = points_start_index + num_points

       contig = modeled_segments_df[s, "CONTIG"]
       offset = contig_starts[match(contig, contig_names)]
       segment_start = offset + modeled_segments_df[s, "START"]
       segment_end = offset + modeled_segments_df[s, "END"]
       genomic_coordinates = offset + denoised_copy_ratios_df[points_start_index:points_end_index, "MIDDLE"]

       denoised_copy_ratios = denoised_copy_ratios_df[points_start_index:points_end_index, "COPY_RATIO"]

       colors = c("coral", "dodgerblue")
       points(x=genomic_coordinates, y=denoised_copy_ratios, col=colors[s %% 2 + 1], pch=".", cex=point_size)

       copy_ratio_posterior_10 = 2^modeled_segments_df[s, "LOG2_COPY_RATIO_POSTERIOR_10"]
       copy_ratio_posterior_50 = 2^modeled_segments_df[s, "LOG2_COPY_RATIO_POSTERIOR_50"]
       copy_ratio_posterior_90 = 2^modeled_segments_df[s, "LOG2_COPY_RATIO_POSTERIOR_90"]
       segments(x0=segment_start, y0=copy_ratio_posterior_50, x1=segment_end, y1=copy_ratio_posterior_50, col="black", lwd=2, lty=1)
       rect(xleft=segment_start, ybottom=copy_ratio_posterior_10, xright=segment_end, ytop=copy_ratio_posterior_90, lwd=1, lty=1)

       points_start_index = points_start_index + num_points
   }
}

PlotAlternateAlleleFractionsWithModeledSegments = function(allelic_counts_df, modeled_segments_df, contig_names, contig_starts, point_size=0.4) {
   points_start_index = 1
   for (s in 1:nrow(modeled_segments_df)) {
       #skip segments with no points
       num_points = modeled_segments_df[s, "NUM_POINTS_ALLELE_FRACTION"]
       if (num_points == 0) {
           next
       }
       points_end_index = points_start_index + num_points

       contig = modeled_segments_df[s, "CONTIG"]
       offset = contig_starts[match(contig, contig_names)]
       segment_start = offset + modeled_segments_df[s, "START"]
       segment_end = offset + modeled_segments_df[s, "END"]
       genomic_coordinates = offset + allelic_counts_df[points_start_index:points_end_index, "POSITION"]

       ref_counts = allelic_counts_df[points_start_index:points_end_index, "REF_COUNT"]
       alt_counts = allelic_counts_df[points_start_index:points_end_index, "ALT_COUNT"]
       alternate_allele_fractions = alt_counts / (alt_counts + ref_counts)

       colors = c("coral", "dodgerblue")
       points(x=genomic_coordinates, y=alternate_allele_fractions, col=colors[s %% 2 + 1], pch=".", cex=point_size)

       minor_allele_fraction_posterior_10 = modeled_segments_df[s, "MINOR_ALLELE_FRACTION_POSTERIOR_10"]
       minor_allele_fraction_posterior_50 = modeled_segments_df[s, "MINOR_ALLELE_FRACTION_POSTERIOR_50"]
       minor_allele_fraction_posterior_90 = modeled_segments_df[s, "MINOR_ALLELE_FRACTION_POSTERIOR_90"]
       segments(x0=segment_start, y0=minor_allele_fraction_posterior_50, x1=segment_end, y1=minor_allele_fraction_posterior_50, col="black", lwd=2, lty=1)
       rect(xleft=segment_start, ybottom=minor_allele_fraction_posterior_10, xright=segment_end, ytop=minor_allele_fraction_posterior_90, lwd=1, lty=1)

       major_allele_fraction_posterior_10 = 1 - minor_allele_fraction_posterior_10
       major_allele_fraction_posterior_50 = 1 - minor_allele_fraction_posterior_50
       major_allele_fraction_posterior_90 = 1 - minor_allele_fraction_posterior_90
       segments(x0=segment_start, y0=major_allele_fraction_posterior_50, x1=segment_end, y1=major_allele_fraction_posterior_50, col="black", lwd=2, lty=1)
       rect(xleft=segment_start, ybottom=major_allele_fraction_posterior_90, xright=segment_end, ytop=major_allele_fraction_posterior_10, lwd=1, lty=1)

       points_start_index = points_start_index + num_points
   }
}
