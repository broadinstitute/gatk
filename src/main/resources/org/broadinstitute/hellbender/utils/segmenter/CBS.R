# Useful for debugging
#options(error = quote({dump.frames(dumpto = "CBS_dump", to.file = TRUE); q()}))

# Library used for segmentation
library(DNAcopy)
library(naturalsort)

library(optparse)
option_list <- list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--targets_file", "-targets_file"), dest="targets_file", action="store"),
    make_option(c("--output_file", "-output_file"), dest="output_file", action="store"),
    make_option(c("--log2_input", "-log"), dest="log2_input", action="store"))

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)
save(opt, file="debug.RData")

sample_name=opt[["sample_name"]]
tn_file=opt[["targets_file"]]
output_file=opt[["output_file"]]
log_input=as.logical(opt[["log2_input"]])

# Use a function for debugging purposes
segment_data = function(sample_name, tn_file, output_file, log_input) {
	# Read in file and extract needed data
	tn = read.table(tn_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	contig = tn[,"contig"]
	pos = tn[,"stop"]
	# Ability to switch between copy-ratio and log2 copy-ratio
	if (log_input) {
	    dat = tn[,sample_name]
	} else {
	    dat = log2(tn[,sample_name])
	}

	# Create CNA object
	cna_dat = CNA(dat, contig, pos, data.type="logratio", sampleid=sample_name)

	# Perform segmentation
	set.seed(25)
	segmented = segment(smooth.CNA(cna_dat))$output

	# Ensure that there are no too-small values which will be problematic for downstream tools.
	segmented[,"seg.mean"] = 2^segmented[,"seg.mean"]
	segmented[segmented[,"seg.mean"]<.Machine$double.eps,"seg.mean"] = .Machine$double.eps

	# Convention for column names
	colnames(segmented) = c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

	# Undo conversion of sample_name to R format (retain dashes, spaces, etc.)
	segmented[,"Sample"] = sample_name

    # Order based on contig (already ordered based on start)
    sorting = unique(naturalsort(segmented[,"Chromosome"]))
    segmented$Chromosome=factor(segmented$Chromosome, levels=sorting)
    segmented = segmented[order(segmented[,"Chromosome"]),]

	# Output seg file
	print(paste("Writing segment file: ", output_file, sep=""))
	write.table(segmented, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
}

segment_data(sample_name, tn_file, output_file, log_input)
