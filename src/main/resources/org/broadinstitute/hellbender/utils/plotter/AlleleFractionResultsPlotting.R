# Useful for debugging
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q()}))

library(optparse)
option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--snp_counts_file", "-snp_counts_file"), dest="snp_counts_file", action="store"),
    make_option(c("--segments_file", "-segments_file"), dest="segments_file", action="store"),
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--sex_chrs", "-sexchrs"), dest="sex_chrs", action="store"))

opt = parse_args(OptionParser(option_list=option_list))
print(opt)
save(opt, file="debug.RData")

sample_name=opt[["sample_name"]]
snp_counts_file=opt[["snp_counts_file"]]
segments_file=opt[["segments_file"]]
output_file=opt[["output_dir"]]
sex_chrs=as.logical(opt[["sex_chrs"]])

## Via http://genome.ucsc.edu/cgi-bin/hgTables
## HG19 chromosome lengths and centromere bounds relative to beginning of chromosomes
chromosome_lengths = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
        146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
        102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566)
centromere_starts = c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679,
	39254935, 51644205, 34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782,
	26369569, 11288129, 13000000)
centromere_ends = c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 50367679,
	42254935, 54644205, 37856694, 19000000, 19000000, 20000000, 38335801, 25263006, 18460898, 27681782,
	29369569, 14288129, 16000000)

###Add X and Y chromosomes###
if(sex_chrs){
    chromosome_lengths = c(chromosome_lengths, 155270560, 59373566)
    centromere_starts = c(centromere_starts, 58632012, 10104553)
    centromere_ends = c(centromere_ends, 61632012, 13104553)
}

#henceforth we use genomic positions scaled relative to the total genome length
genome_length = sum(chromosome_lengths)
chromosome_lengths = chromosome_lengths / genome_length
centromere_starts = centromere_starts / genome_length   #relative to chromosome starts
centromere_ends = centromere_ends / genome_length       #relative to chromosome ends

num_chromosomes = length(chromosome_lengths)

chromosome_ends = cumsum(chromosome_lengths)
chromosome_starts = c(0, head(chromosome_ends, -1))
chromosome_centers = (chromosome_starts + chromosome_ends) / 2

centromere_starts = chromosome_starts + centromere_starts
centromere_ends = chromosome_starts + centromere_ends

# Use a function for debugging purposes
create_allele_fraction_plots_file = function(sample_name, snp_counts_file, segments_file, output_dir) {
	# Read in file and extract needed data
	snp_counts = read.table(snp_counts_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	segments = read.table(segments_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)


	snpContig = snp_counts[,"CONTIG"]
    snpContig[snpContig=="X"]=23
    snpContig[snpContig=="Y"]=24
    snp_counts[,"CONTIG"]=as.numeric(snpContig)

    segments[segments["Chromosome"]=="X","Chromosome"]=23
    segments[segments["Chromosome"]=="Y","Chromosome"]=24
    segments[,"Chromosome"]=as.numeric(segments[,"Chromosome"])

    #exclude segments with no SNPS
    segments = segments[segments$Num_SNPs > 0,]

    plot_file_name = file.path( output_dir, paste(sample_name, "_FullGenome", ".png", sep=""))
    png(plot_file_name, 12, 7, units="in", type="cairo", res=300, bg="white")
    SetUpPlot("Allele Fraction", 0, 1, "Chromosome", TRUE)
    AlleleFractionPlot(snp_counts, segments)
    dev.off()
}

## Chromosome Background for data
SetUpPlot = function(y.lab, y.min, y.max, x.lab, label_chromosomes){
  suppressWarnings( par( mar=c(3.1,3.6,0.1,0),mgp=c(2,-0.2,-1.1) ) )
  plot(0, type = "n", bty = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "", ylab = "", main = "", xaxt = "n")
  mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
  mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)

  if(label_chromosomes){
  	labels = c(1:num_chromosomes)
  	odd_indices = labels %% 2 == 1
  	mtext(text = labels[odd_indices], side = 1, line = -0.45, at = chromosome_centers[odd_indices],
    		las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  	mtext(text = labels[!odd_indices], side = 1, line = 0, at = chromosome_centers[!odd_indices],
      	las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  }

  for (i in 1:(num_chromosomes - 1)) {
    use.col = ifelse(i%%2 == 1, "grey90", "white")
    rect(xleft = chromosome_starts[i], ybottom = y.min, xright = chromosome_ends[i],
    		ytop = y.max, col = use.col, border = NA)
    lines(y = c(y.min, y.max), x = rep(centromere_starts[i], 2), lty = 3, lwd = 0.5)
    lines(y = c(y.min, y.max), x = rep(centromere_ends[i], 2), lty = 3, lwd = 0.5)
  }
}

AlleleFractionPlot = function(snp_df, segments_df){
  for( s in 1:nrow(segments_df)) {
     snp_indices=which(snp_df[,"CONTIG"]==segments_df[s,"Chromosome"] & snp_df[,"POS"] > segments_df[s,"Start"] & snp_df[,"POS"] < segments_df[s,"End"])
     genomic_coordinates = chromosome_starts[snp_df[snp_indices,"CONTIG"]] + snp_df[snp_indices,"POS"] / genome_length
     ref_counts = snp_df[snp_indices, "REF_COUNT"]
     alt_counts = snp_df[snp_indices, "ALT_COUNT"]
     MAF = ifelse(alt_counts<ref_counts, alt_counts/(alt_counts+ref_counts), ref_counts/(alt_counts+ref_counts))

     colors=c("coral", "dodgerblue")
     points(x=genomic_coordinates, y=MAF, col=colors[s%%2+1], pch=16, cex=0.2)

     lines(x=range(genomic_coordinates), y=rep(segments_df[s, "MAF_Post_Mean"], 2), col="black", lwd=1, lty=2)

     #draw error bars
     segment_mean = segments_df[s, "MAF_Post_Mean"]
     segment_std = segments_df[s, "MAF_Post_Std"]
     middle_of_segment = (max(genomic_coordinates) + min(genomic_coordinates))/2
     arrows(middle_of_segment,segment_mean - segment_std,middle_of_segment,segment_mean + segment_std, code=3, length=0.02, angle = 90)
  }
}

create_allele_fraction_plots_file(sample_name, snp_counts_file, segments_file, output_file)


