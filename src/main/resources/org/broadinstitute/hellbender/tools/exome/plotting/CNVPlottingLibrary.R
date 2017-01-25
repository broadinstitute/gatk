## Contig Background for data
SetUpPlot = function(y.lab, y.min, y.max, x.lab, contig_names, contig_starts, contig_ends, do_label_contigs) {
    num_contigs = length(contig_names)
    contig_centers = (contig_starts + contig_ends) / 2
    genome_length = contig_ends[num_contigs]
    suppressWarnings(par(mar=c(3.1, 3.6, 0.1, 0), mgp=c(2, -0.2, -1.1)))
    plot(0, type="n", bty="n", xlim=c(0, genome_length), ylim=c(0, y.max), xlab="", ylab="", main="", xaxt="n")
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

PlotCopyRatio = function(coverage_df, color, contig_names, contig_starts) {
    genomic_coordinates = contig_starts[match(coverage_df$contig, contig_names)] + coverage_df$stop
    points(x=genomic_coordinates, y=coverage_df$VALUE, col=color, pch=16, cex=0.2)
}

PlotCopyRatioWithSegments = function(coverage_df, segments_df, contig_names, contig_starts, is_ACNV, point_size=0.2) {
    for (s in 1:nrow(segments_df)) {
         contig = segments_df[s, "Chromosome"]
         segment_start = segments_df[s, "Start"]
         segment_end = segments_df[s, "End"]

         target_indices = which(coverage_df[, "contig"] == contig & coverage_df[, "stop"] > segment_start & coverage_df[, "stop"] < segment_end)

         offset = contig_starts[match(contig, contig_names)]
         genomic_coordinates = offset + coverage_df[target_indices, "stop"]

         colors = c("coral", "dodgerblue")
         points(x=genomic_coordinates, y=coverage_df[target_indices, "VALUE"], col=colors[s %% 2 + 1], pch=16, cex=point_size)

         if (is_ACNV) {
             segment_mode = 2^segments_df[s, "Segment_Mean_Post_Mode"]
             segment_high = 2^segments_df[s, "Segment_Mean_Post_Hi"]
             segment_low = 2^segments_df[s, "Segment_Mean_Post_Lo"]
             segments(x0=min(genomic_coordinates), y0=segment_mode, x1=max(genomic_coordinates), y1=segment_mode, col="black", lwd=2, lty=1)
             rect(xleft=min(genomic_coordinates), ybottom=segment_low, xright=max(genomic_coordinates), ytop=segment_high, lwd=1, lty=1)
         } else {
             segment_mean = segments_df[s, "Segment_Mean"]
             segments(x0=min(genomic_coordinates), y0=segment_mean, x1=max(genomic_coordinates), y1=segment_mean, col="black", lwd=2, lty=1)
        }
    }
}

PlotAlleleFractionWithSegments = function(snp_df, segments_df, contig_names, contig_starts, point_size=0.2) {
    for (s in 1:nrow(segments_df)) {
         #skip segments with no SNPs
         if (segments_df[s, "Num_SNPs"] == 0) {
            next
         }
         contig = segments_df[s, "Chromosome"]
         segment_start = segments_df[s, "Start"]
         segment_end = segments_df[s, "End"]
         snp_indices = which(snp_df[, "CONTIG"] == contig & snp_df[,"POSITION"] > segment_start & snp_df[, "POSITION"] < segment_end)

         offset = contig_starts[match(contig, contig_names)]
         genomic_coordinates = offset + snp_df[snp_indices, "POSITION"]

         ref_counts = snp_df[snp_indices, "REF_COUNT"]
         alt_counts = snp_df[snp_indices, "ALT_COUNT"]
         MAF = ifelse(alt_counts < ref_counts, alt_counts / (alt_counts + ref_counts), ref_counts / (alt_counts + ref_counts))

         colors = c("coral", "dodgerblue")
         points(x=genomic_coordinates, y=MAF, col=colors[s %% 2 + 1], pch=16, cex=point_size)

         segment_mode = segments_df[s, "MAF_Post_Mode"]
         segment_high = segments_df[s, "MAF_Post_Hi"]
         segment_low = segments_df[s, "MAF_Post_Lo"]
         segments(x0=min(genomic_coordinates), y0=segment_mode, x1=max(genomic_coordinates), y1=segment_mode, col="black", lwd=2, lty=1)
         rect(xleft=min(genomic_coordinates), ybottom=segment_low, xright=max(genomic_coordinates), ytop=segment_high, lwd=1, lty=1)
    }
}