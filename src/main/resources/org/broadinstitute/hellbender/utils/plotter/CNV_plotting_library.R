## Via http://genome.ucsc.edu/cgi-bin/hgTables
## HG19 chromosome lengths (autosomal 1-22, followed by X and Y) and centromere bounds relative to beginning of chromosomes
chromosome_lengths = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
        146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
        102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
centromere_starts = c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679,
	39254935, 51644205, 34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782,
	26369569, 11288129, 13000000, 58632012, 10104553)
centromere_ends = c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 50367679,
	42254935, 54644205, 37856694, 19000000, 19000000, 20000000, 38335801, 25263006, 18460898, 27681782,
	29369569, 14288129, 16000000, 61632012, 13104553)

#henceforth we use genomic positions scaled relative to the total genome length

chromosome_ends = cumsum(chromosome_lengths)
chromosome_starts = c(0, head(chromosome_ends, -1))
chromosome_centers = (chromosome_starts + chromosome_ends) / 2
centromere_starts = chromosome_starts + centromere_starts
centromere_ends = chromosome_starts + centromere_ends
chromosome_labels = c(1:22, "X","Y")

convert_XY_to_23_24 = function(data) {
    data[data=="X"] = 23
    data[data=="Y"] = 24
    as.numeric(data)
}

## Chromosome Background for data
SetUpPlot = function(y.lab, y.min, y.max, x.lab, label_chromosomes, num_chromosomes){
  genome_length = sum(chromosome_lengths[1:num_chromosomes])
  suppressWarnings( par( mar=c(3.1,3.6,0.1,0),mgp=c(2,-0.2,-1.1) ) )
  plot(0, type = "n", bty = "n", xlim = c(0, genome_length), ylim = c(0, y.max), xlab = "", ylab = "", main = "", xaxt = "n")
  mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
  mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)

  if(label_chromosomes){
  	mtext(text = chromosome_labels[1:num_chromosomes], side = 1, line = ifelse(c(1:num_chromosomes)%%2 == 1, -0.45, 0.0),
  	at = chromosome_centers[1:num_chromosomes], las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  }

  for (i in 1:num_chromosomes) {
    use.col = ifelse(i%%2 == 1, "grey90", "white")
    rect(xleft = chromosome_starts[i], ybottom = y.min, xright = chromosome_ends[i],
    		ytop = y.max, col = use.col, border = NA)
    lines(y = c(y.min, y.max), x = rep(centromere_starts[i], 2), lty = 3, lwd = 0.5)
    lines(y = c(y.min, y.max), x = rep(centromere_ends[i], 2), lty = 3, lwd = 0.5)
  }
}

#this takes in SNPs as input to ACNV and segments as output by ACNV
PlotAlleleFraction = function(snp_df, segments_df){
  for( s in 1:nrow(segments_df)) {
     #skip segments with no SNPs
     if (segments_df[s, "Num_SNPs"] == 0) {
        next
     }
     contig = segments_df[s,"Chromosome"]
     segment_start = segments_df[s, "Start"]
     segment_end = segments_df[s, "End"]
     snp_indices=which(snp_df[,"CONTIG"] == contig & snp_df[,"POSITION"] > segment_start & snp_df[,"POSITION"] < segment_end)
     genomic_coordinates = chromosome_starts[contig] + snp_df[ snp_indices,"POSITION"]
     ref_counts = snp_df[snp_indices, "REF_COUNT"]
     alt_counts = snp_df[snp_indices, "ALT_COUNT"]
     MAF = ifelse(alt_counts < ref_counts, alt_counts/(alt_counts + ref_counts), ref_counts/(alt_counts + ref_counts))

     colors=c("coral", "dodgerblue")
     points(x=genomic_coordinates, y=MAF, col=colors[s%%2+1], pch=16, cex=0.2)

     segment_mode = segments_df[s, "MAF_Post_Mode"]
     segment_high = segments_df[s, "MAF_Post_Hi"]
     segment_low = segments_df[s, "MAF_Post_Lo"]
     segments(x0=min(genomic_coordinates), y0=segment_mode, x1=max(genomic_coordinates), y1=segment_mode, col="black", lwd=2, lty=1)
     rect(xleft=min(genomic_coordinates), ybottom=segment_low, xright=max(genomic_coordinates), ytop=segment_high, lwd=1, lty=1)
  }
}

PlotCopyRatio = function(df, color) {
  genome.crds = chromosome_starts[df$contig] + df$stop
  points(x=genome.crds, y=df$VALUE, col=color, pch=16, cex=0.2)
}

PlotCopyRatioWithSegments = function(coverage_df, segments_df, is_ACNV ) {
  for(s in 1:nrow(segments_df)) {
     contig = segments_df[s,"Chromosome"]
     segment_start = segments_df[s, "Start"]
     segment_end = segments_df[s, "End"]

     target_indices=which(coverage_df[,"contig"] == contig & coverage_df[,"stop"] > segment_start & coverage_df[,"stop"] < segment_end)
     genomic_coordinates = chromosome_starts[contig] + coverage_df[ target_indices,"stop"]
     colors=c("coral", "dodgerblue")
     points(x=genomic_coordinates, y=coverage_df[target_indices, "VALUE"], col=colors[s%%2+1], pch=16, cex=0.2)

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