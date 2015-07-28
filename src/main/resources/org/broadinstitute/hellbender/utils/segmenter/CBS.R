# Useful for debugging
#options(error = quote({dump.frames(dumpto = "CBS_dump", to.file = TRUE); q()}))

# Library used for segmentation
library(DNAcopy)
library(naturalsort)

# Parse arguments
library(argparser)
p <- arg.parser("A text file modifying program");
p <- add.argument(p, "sample_name", help="sample name");
p <- add.argument(p, "targets_file", help="input file");
p <- add.argument(p, "output_file", help="output file name");
p <- add.argument(p, "min_log_value", help="values will be thresholded to this value", default=-10);
argv = parse.args(p, argv = commandArgs(trailingOnly = TRUE))
print(argv)

sample_name=argv[["sample_name"]]
tn_file=argv[["targets_file"]]
output_file=argv[["output_file"]]

# Use a function for debugging purposes
segment_data = function(sample_name, tn_file, output_file, min_log_value=-10) {
	# Read in file and extract needed data
	tn = read.table(tn_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	contig = tn[,"CONTIG"]
	pos = tn[,"END"]
	dat = 2^tn[,sample_name]

	# Create CNA object
	cna_dat = CNA(dat, contig, pos, data.type="logratio", sampleid=sample_name)

	# Perform segmentation
	set.seed(25)
	segmented = segment(smooth.CNA(cna_dat))$output
	segmented[,"seg.mean"] = log2(segmented[,"seg.mean"])

	# Ensure that there are no too-small values which will be problematic for downstream tools.
	log10_min_copy_ratio = min_log_value
	segmented[segmented[,"seg.mean"]<log10_min_copy_ratio,"seg.mean"] = log10_min_copy_ratio

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

segment_data(sample_name, tn_file, output_file)