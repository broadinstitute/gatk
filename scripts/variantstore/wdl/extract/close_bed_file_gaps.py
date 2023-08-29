import csv
import argparse
import re
import math

# put these up here to get the defaults.  But override them later with user input
binsize_kb = 100
infile = f"gvs_vet_weights_{binsize_kb}kb.bed"
outfile = f"gvs_vet_weights_{binsize_kb}kb_gaps_filled.bed"

argParser = argparse.ArgumentParser()

argParser.add_argument("-s", "--bin-size", help="size of the bins in kb", type=int, default=binsize_kb)
argParser.add_argument("-i", "--infile", help="bed input file", type=str, default=infile)
argParser.add_argument("-o", "--outfile", help="bed output file", type=str, default=infile)
argParser.add_argument("-l", "--interval-list", help="an interval list we can use to fill in any data on a contig past the end of our samples", type=str)

args = argParser.parse_args()

binsize_kb = args.bin_size
infile = args.infile
outfile = args.outfile

contigSizes = {}
# an interval list is not required.  But if one is supplied, we'll use it to determine how long each contig is so we can
# pad out our data to that size.  We yank the contig name and length info directly from the lines using simple regexes
if args.interval_list:
	interval_list_file = open(args.interval_list, 'r')
	for header_line in interval_list_file:
		if header_line.startswith("@SQ"):
			# look for the contig and length fields
			contig_result = re.search(r'SN:(\S+)', header_line);
			length_result = re.search(r'LN:(\d+)', header_line);
		
			if contig_result and length_result:
				contig = contig_result.group(1)
				contig_length = length_result.group(1)
				contigSizes[contig] = int(contig_length) 			
else:
	print("We do not have an interval list")



lastContig = ""
lastEnding = 0
intervalSize = binsize_kb * 1000
print(f"Interval size is {intervalSize}")

with open(infile, newline='') as bedfile, open(outfile, 'w') as output:
	bedreader = csv.reader(bedfile, delimiter='\t')
	# logic is simple.  For every row, print it to the new file.  EXCEPT, we also want to 
	# detect gaps in between intervals and fill them with 0 weight intervals.  That way, 
	# when it's read in we'll have a solid block of data for an entire contig and can put
	# it in an array for fast looking
	for row in bedreader:
		contig = row[0]
		start = int(row[1])
		ending = int(row[2])
		# Logic to detect gaps in the intervals is only valid WITHIN a contig
		if contig == lastContig:
			# the intervals are read as [,) so the final entry in the previous one should
			# match the first entry in the next one if there are no gaps
			if start != lastEnding:
				print(f"detected gap from {contig} {lastEnding} to {start}")
				sizeOfGap = start - lastEnding
				numEmptyEntries = sizeOfGap / intervalSize
				print(f"Will create {numEmptyEntries} 0 weight entries here")
				for s in range(lastEnding, start, intervalSize):
					print(f"\t{contig} {s} - {s + intervalSize}")
					newRow = [contig, str(s), str(s+intervalSize), ".", "0"]
					output.write("\t".join(newRow) + "\n")
		else:
			# see if there's anything on the end of the last contig to extend out to
			if lastContig in contigSizes:
				# we'll want to round up the listed ending to the next block in case we
				# haven't reached it yet
				newEndingBlock = int(math.ceil(contigSizes[lastContig] / float(intervalSize))) * intervalSize
				print(f"Ending detected for contig: {lastContig} rounded up to block {newEndingBlock}.  Last block written was {lastEnding}")
				for s in range(lastEnding, newEndingBlock, intervalSize):
					print(f"\t{lastContig} {s} - {s + intervalSize}")
					newRow = [lastContig, str(s), str(s+intervalSize), ".", "0"]
					output.write("\t".join(newRow) + "\n")
			# make sure we always start at 0 on a new contig and search for gaps
			print(f"detected gap at beginning of {contig}")
			lastEnding = 0
			sizeOfGap = start
			numEmptyEntries = sizeOfGap / intervalSize
			print(f"Would create {numEmptyEntries} 0 weight entries here")
			for s in range(lastEnding, start, intervalSize):
				print(f"\t{contig} {s} - {s + intervalSize}")
				newRow = [contig, str(s), str(s+intervalSize), ".", "0"]
				output.write("\t".join(newRow) + "\n")
			
		lastContig = contig
		lastEnding = ending
		# now that we've inserted any fake rows necessary to fill gaps, write the current one
		output.write("\t".join(row) + "\n")
		
		
	#end for loop
	# check for anything leftover at the end of the last contig
	if lastContig in contigSizes:
		# we'll want to round up the listed ending to the next block in case we
		# haven't reached it yet
		newEndingBlock = int(math.ceil(contigSizes[lastContig] / float(intervalSize))) * intervalSize
		print(f"Ending detected for contig: {lastContig} rounded up to block {newEndingBlock}.  Last block written was {lastEnding}")
		for s in range(lastEnding, newEndingBlock, intervalSize):
			print(f"\t{lastContig} {s} - {s + intervalSize}")
			newRow = [lastContig, str(s), str(s+intervalSize), ".", "0"]
			output.write("\t".join(newRow) + "\n")
	print("Processing done")


