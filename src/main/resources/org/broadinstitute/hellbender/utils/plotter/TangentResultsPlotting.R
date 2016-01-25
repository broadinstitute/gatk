# Useful for debugging
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q()}))

library(optparse)

option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--targets_file", "-targets_file"), dest="targets_file", action="store"),
    make_option(c("--pre_tn_file", "-pre_tn_file"), dest="pre_tn_file", action="store"),
    make_option(c("--segments_file", "-segments_file"), dest="seg_file", action="store"),
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--log2_input", "-log"), dest="log2_input", action="store"),
    make_option(c("--sex_chrs", "-sexchrs"), dest="sex_chrs", action="store"))

opt = parse_args(OptionParser(option_list=option_list))
print(opt)
save(opt, file="debug.RData")

sample_name=opt[["sample_name"]]
tn_file=opt[["targets_file"]]
pre_tn_file=opt[["pre_tn_file"]]
segments_file=opt[["seg_file"]]
output_file=opt[["output_dir"]]
log_input=as.logical(opt[["log2_input"]])
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
create_tangent_plots_file = function(sample_name, tn_file, pre_tn_file, segments_file, output_dir, log_input, sex_chrs) {
	# Read in file and extract needed data
	segments = read.table(segments_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	tn = read.table(tn_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	preTn = read.table(pre_tn_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)

	# Ability to switch between copy-ratio and log2 copy-ratio
	if (log_input) {
	    tnDat = 2^tn[,sample_name]
	    preTnDat = 2^preTn[,sample_name]
	}else {
	    tnDat = tn[,sample_name]
	    preTnDat = preTn[,sample_name]
	}

    tnContig = tn[,"contig"]
    tnContig[tnContig=="X"]=23
    tnContig[tnContig=="Y"]=24
    tnContig=as.numeric(tnContig)
    tnPos = tn[,"stop"]

    preTnContig = preTn[,"contig"]
    preTnContig[preTnContig=="X"]=23
    preTnContig[preTnContig=="Y"]=24
    preTnContig=as.numeric(preTnContig)
    preTnPos = preTn[,"stop"]

    segments[segments["Chromosome"]=="X","Chromosome"]=23
    segments[segments["Chromosome"]=="Y","Chromosome"]=24
    segments[,"Chromosome"]=as.numeric(segments[,"Chromosome"])

    preQc = QC(preTnDat)
    postQc = QC(tnDat)

    write.table(round(preQc, 3), file.path( output_dir, paste(sample_name, "_preQc.txt", sep="")),
    	quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(round(postQc, 3), file.path( output_dir, paste(sample_name, "_postQc.txt", sep="")),
    	quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(round(preQc-postQc,3), file.path( output_dir, paste(sample_name, "_dQc.txt", sep="")),
    	quote=FALSE, row.names=FALSE, col.names=FALSE)

    plot.fn = file.path( output_dir, paste(sample_name, "_FullGenome", ".png", sep=""))
    png(plot.fn, 12, 7, units="in", type="cairo", res=300, bg="white")
    SetUpPlot("Total copy ratio", 0, 4, "", TRUE, sex_chrs)
    TotalCopyRatioWithSegments(tnDat, tnPos, tnContig, segments)
    dev.off()

    pre_color_blue="#3B5DFF"
    post_color_green="#4FC601"
    plot.fn = file.path( output_dir, paste(sample_name, "_Before_After_CR_Lim_4", ".png", sep=""))
    png(plot.fn, 12, 7, units="in", type="cairo", res=300, bg="white")
    par( mfrow=c(2,1), cex=0.75, las=1 )
    SetUpPlot("Copy-Ratio", 0, 4, paste("Pre Tangent Normalization, pre-Qc = ", round(preQc, 3), sep=""), FALSE, sex_chrs)
    TotalCopyRatio(preTnDat, preTnPos, preTnContig, sample_name, pre_color_blue)
    SetUpPlot("Copy-Ratio", 0, 4, paste("Post Tangent Normalization, post-Qc = ", round(postQc, 3), ", dQc = ",
        round(preQc-postQc,3), ", scaled_Qc = ", round((preQc-postQc)/preQc,3), sep=""), TRUE, sex_chrs)
    TotalCopyRatio(tnDat, tnPos, tnContig, sample_name, post_color_green)
    dev.off()

    plot.fn = file.path( output_dir, paste(sample_name, "_Before_After", ".png", sep=""))
    png(plot.fn, 12, 7, units="in", type="cairo", res=300, bg="white")
    par( mfrow=c(2,1), cex=0.75, las=1 )
    SetUpPlot("Copy-Ratio", 0, max(preTnDat), paste("Pre Tangent Normalization, pre-Qc = ", round(preQc, 3), sep=""), FALSE, sex_chrs)
    TotalCopyRatio(preTnDat, preTnPos, preTnContig, sample_name, pre_color_blue)
    SetUpPlot("Copy-Ratio", 0, max(tnDat), paste("Post Tangent Normalization, post-Qc = ", round(postQc, 3),
        ", dQc = ", round(preQc-postQc,3), ", scaled_dQc = ", round((preQc-postQc)/preQc,3), sep=""), TRUE, sex_chrs)
    TotalCopyRatio(tnDat, tnPos, tnContig, sample_name, post_color_green)
    dev.off()
}

QC = function(dat) {
    return(median(abs(diff(dat))))
}

## Chromosome Background for data
SetUpPlot = function(y.lab, y.min, y.max, x.lab, label_chromosomes, sex_chrs){
  suppressWarnings( par( mar=c(3.1,3.6,0.1,0),mgp=c(2,-0.2,-1.1) ) )
  plot(0, type = "n", bty = "n", xlim = c(0, 1), ylim = c(y.min, y.max),
       xlab = "", ylab = "", main = "", xaxt = "n")
  mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
  mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)

  if(label_chromosomes){
  	labels = c(1:num_chromosomes)
  	odd_indices = labels %% 2 == 1
  	if(sex_chrs){
        labels[23]="X"
        labels[24]="Y"
    }

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

TotalCopyRatioWithSegments = function(dat, pos, contig, segments){
  for( s in 1:nrow(segments)) {
     segContig=segments[s,"Chromosome"]
     End=segments[s,"End"]
     if(s>1){
        last_End=segments[s-1,"End"]
     }
     if(s==1||segContig!=segments[s-1,"Chromosome"]){
        last_End=segments[s,"Start"]
     }
     seg.mean = segments[s,"Segment_Mean"]
     segment_indices = which(contig==segContig & pos>last_End & pos<=End)
     cur_pos=pos[segment_indices]
     cur_dat=dat[segment_indices]
     genome.crds = chromosome_starts[segContig] + cur_pos / genome_length
     colors=c("coral", "dodgerblue")
     points(genome.crds, cur_dat, col=colors[s%%2+1], pch=16, cex=0.2)
     lines(range(genome.crds), rep(seg.mean, 2), col="black", lwd=1, lty=2)
  }
}

TotalCopyRatio = function(dat, pos, contig, sample_name, color){
  genome.crds = chromosome_starts[contig] + pos / genome_length
  points(genome.crds, dat, col=color, pch=16, cex=0.2)
}

create_tangent_plots_file(sample_name, tn_file, pre_tn_file, segments_file, output_file, log_input, sex_chrs)