# Useful for debugging
#options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q()}))

library(optparse)

option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--targets_file", "-targets_file"), dest="targets_file", action="store"),
    make_option(c("--pre_tn_file", "-pre_tn_file"), dest="pre_tn_file", action="store"),
    make_option(c("--seg_file", "-seg_file"), dest="seg_file", action="store"),
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--log2_input", "-log"), dest="log2_input", action="store"))

opt = parse_args(OptionParser(option_list=option_list))
print(opt)
save(opt, file="debug.RData")

sample_name=opt[["sample_name"]]
tn_file=opt[["targets_file"]]
pre_tn_file=opt[["pre_tn_file"]]
seg_file=opt[["seg_file"]]
output_file=opt[["output_dir"]]
log_input=as.logical(opt[["log2_input"]])

# Use a function for debugging purposes
create_tangent_plots_file = function(sample_name, tn_file, pre_tn_file, seg_file, output_dir, log_input) {
	# Read in file and extract needed data
	seg = read.table(seg_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
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
    tnPos = tn[,"stop"]
    preTnContig = preTn[,"contig"]
    preTnPos = preTn[,"stop"]
    preQc = QC(preTnDat)
    postQc = QC(tnDat)

    plot.fn = file.path( output_dir, paste(sample_name, "_FullGenome", ".png", sep=""))
    png(plot.fn, 12, 7, units="in", type="cairo", res=300, bg="white")
    chr.dat = SetUpPlot("Total copy ratio", 0, 4, "", TRUE)
    TotalCopyRatioWithSegments(tnDat, tnPos, tnContig, seg, chr.dat)
    dev.off()

    pre_color_blue="#3B5DFF"
    post_color_green="#4FC601"
    plot.fn = file.path( output_dir, paste(sample_name, "_Before_After_CR_Lim_4", ".png", sep=""))
    png(plot.fn, 12, 7, units="in", type="cairo", res=300, bg="white")
    par( mfrow=c(2,1), cex=0.75, las=1 )
    SetUpPlot("Copy-Ratio", 0, 4, paste("Pre Tangent Normalization, pre-Qc = ", round(preQc, 3), sep=""), FALSE)
    TotalCopyRatio(preTnDat, preTnPos, preTnContig, chr.dat, sample_name, pre_color_blue)
    SetUpPlot("Copy-Ratio", 0, 4, paste("Post Tangent Normalization, post-Qc = ", round(postQc, 3), ", dQc = ",
        round(preQc-postQc,3), ", scaled_Qc = ", round((preQc-postQc)/preQc,3), sep=""), TRUE)
    TotalCopyRatio(tnDat, tnPos, tnContig, chr.dat, sample_name, post_color_green)
    dev.off()

    plot.fn = file.path( output_dir, paste(sample_name, "_Before_After", ".png", sep=""))
    png(plot.fn, 12, 7, units="in", type="cairo", res=300, bg="white")
    par( mfrow=c(2,1), cex=0.75, las=1 )
    SetUpPlot("Copy-Ratio", 0, max(preTnDat), paste("Pre Tangent Normalization, pre-Qc = ", round(preQc, 3), sep=""), FALSE)
    TotalCopyRatio(preTnDat, preTnPos, preTnContig, chr.dat, sample_name, pre_color_blue)
    SetUpPlot("Copy-Ratio", 0, max(tnDat), paste("Post Tangent Normalization, post-Qc = ", round(postQc, 3),
        ", dQc = ", round(preQc-postQc,3), ", scaled_dQc = ", round((preQc-postQc)/preQc,3), sep=""), TRUE)
    TotalCopyRatio(tnDat, tnPos, tnContig, chr.dat, sample_name, post_color_green)
    dev.off()
}

QC = function(dat) {
    return(median(abs(diff(dat))))
}

## Chromosome Background for data
SetUpPlot = function(y.lab, y.min, y.max, x.lab, lab.chr){
  ## Via http://genome.ucsc.edu/cgi-bin/hgTables
  ## HG19 chromosome lengths
  chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
        146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
        102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566)

  ## HG19 centromere positions left and right bounds
  cent.posS = c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679,
	39254935, 51644205, 34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782,
	26369569, 11288129, 13000000)
  cent.posE = c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 50367679,
	42254935, 54644205, 37856694, 19000000, 19000000, 20000000, 38335801, 25263006, 18460898, 27681782,
	29369569, 14288129, 16000000)

  ###Add sex chromosomes###
  ## Add X - DISABLED
  #chr.lens = c(chr.lens, 155270560)
  #cent.posS = c(cent.posS, 58632012)
  #cent.posE = c(cent.posE, 61632012)

  ## Add Y - DISABLED
  #chr.lens = c(chr.lens, 59373566)
  #cent.posS = c(cent.posS, 10104553)
  #cent.posE = c(cent.posE, 13104553)

  genome.len = sum(chr.lens)
  chr.w = chr.lens / genome.len
  suppressWarnings( par( mar=c(3.1,3.6,0.1,0),mgp=c(2,-0.2,-1.1) ) )
  plot(0, type = "n", bty = "n", xlim = c(0, 1), ylim = c(y.min, y.max),
       xlab = "", ylab = "", main = "", xaxt = "n")
  mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
  mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)
  chr.mids=cumsum(chr.w)-chr.w/2
  if(lab.chr){
  	lab.vals = (c(1:length(chr.w)))
  	odd.ix = lab.vals %% 2 == 1
  	mtext(text = lab.vals[odd.ix], side = 1, line = -0.45, at = chr.mids[odd.ix],
    		las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  	mtext(text = lab.vals[!odd.ix], side = 1, line = 0, at = chr.mids[!odd.ix],
      	las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  }
  chr.offsets = c(0, cumsum(chr.w[c(1:(length(chr.w) - 1))]), 1)
  gen.posS = cent.posS/genome.len + chr.offsets[1:(length(chr.offsets) - 1)]
  gen.posE = cent.posE/genome.len + chr.offsets[1:(length(chr.offsets) - 1)]
  for (i in 1:(length(chr.offsets) - 1)) {
    use.col = ifelse(i%%2 == 1, "grey90", "white")
    rect(xleft = chr.offsets[i], ybottom = y.min, xright = chr.offsets[i + 1],
    		ytop = y.max, col = use.col, border = NA)
    lines(y = c(y.min, y.max), x = rep(gen.posS[i], 2), lty = 3, lwd = 0.5)
    lines(y = c(y.min, y.max), x = rep(gen.posE[i], 2), lty = 3, lwd = 0.5)
  }
  chr.dat = list(chr.w, chr.lens, chr.offsets)
  return( chr.dat )
}

TotalCopyRatioWithSegments = function(dat, pos, contig, seg, chr.dat){
  chr.w = chr.dat[[1]]
  chr.lens = chr.dat[[2]]
  chr.offsets = chr.dat[[3]]
  for( s in 1:nrow(seg)) {
     segContig=seg[s,"Chromosome"]
     End=seg[s,"End"]
     if(s>1){
        last_End=seg[s-1,"End"]
     }
     if(s==1||segContig!=seg[s-1,"Chromosome"]){
        last_End=seg[s,"Start"]
     }
     seg.mean=seg[s,"Segment_Mean"]
     seg_ind=which(contig==segContig&pos>last_End&pos<=End)
     cur_pos=pos[seg_ind]
     cur_dat=dat[seg_ind]
     genome.crds = chr.offsets[segContig] + cur_pos / chr.lens[segContig] * chr.w[segContig]
     colors=c("coral", "dodgerblue")
     points(genome.crds, cur_dat, col=colors[s%%2+1], pch=16, cex=0.2)
     lines(range(genome.crds), rep(seg.mean, 2), col="black", lwd=1, lty=2)
  }
}

TotalCopyRatio = function(dat, pos, contig, chr.dat, sample_name, color){
  chr.w = chr.dat[[1]]
  chr.lens = chr.dat[[2]]
  chr.offsets = chr.dat[[3]]
  genome.crds = chr.offsets[contig] + pos / chr.lens[contig] * chr.w[contig]
  points(genome.crds, dat, col=color, pch=16, cex=0.2)
}

create_tangent_plots_file(sample_name, tn_file, pre_tn_file, seg_file, output_file, log_input)
