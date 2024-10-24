import csv
import argparse
import re
import math

argParser = argparse.ArgumentParser()

argParser.add_argument("-s", "--bin-size", help="size of the bins in kb", type=int)
argParser.add_argument("-i", "--infile", help="bed input file", type=str)
argParser.add_argument("-o", "--outfile", help="bed output file", type=str)
argParser.add_argument("-l", "--interval-list", help="an interval list we can use to fill in any data on a contig past the end of our samples", type=str)

args = argParser.parse_args()

binsize_kb = args.bin_size
infile = args.infile
outfile = args.outfile

contig_sizes = {}
# an interval list is not required.  But if one is supplied, we'll use it to determine how long each contig is so we can
# pad out our data to that size.  We yank the contig name and length info directly from the lines using simple regexes
if args.interval_list:
	with open(args.interval_list, 'r') as interval_list_file:
		for header_line in interval_list_file:
			if header_line.startswith("@SQ"):
				# look for the contig and length fields
				contig_result = re.search(r'SN:(\S+)', header_line)
				length_result = re.search(r'LN:(\d+)', header_line)

				if contig_result and length_result:
					contig = contig_result.group(1)
					contig_length = length_result.group(1)
					contig_sizes[contig] = int(contig_length)
else:
	print("We do not have an interval list")

last_contig = ""
last_ending = 0
interval_size = binsize_kb * 1000
print(f"Interval size is {interval_size}")

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
		if contig == last_contig:
			# the intervals are read as [,) so the final entry in the previous one should
			# match the first entry in the next one if there are no gaps
			if start != last_ending:
				print(f"detected gap from {contig} {last_ending} to {start}")
				size_of_gap = start - last_ending
				num_empty_entries = size_of_gap / interval_size
				print(f"Will create {num_empty_entries} 0 weight entries here")
				for s in range(last_ending, start, interval_size):
					print(f"\t{contig} {s} - {s + interval_size}")
					new_row = [contig, str(s), str(s + interval_size), ".", "0"]
					output.write("\t".join(new_row) + "\n")
		else:
			# see if there's anything on the end of the last contig to extend out to
			if last_contig in contig_sizes:
				# we'll want to round up the listed ending to the next block in case we
				# haven't reached it yet
				new_ending_block = int(math.ceil(contig_sizes[last_contig] / float(interval_size))) * interval_size
				print(f"Ending detected for contig: {last_contig} rounded up to block {new_ending_block}.  Last block written was {last_ending}")
				for s in range(last_ending, new_ending_block, interval_size):
					print(f"\t{last_contig} {s} - {s + interval_size}")
					new_row = [last_contig, str(s), str(s + interval_size), ".", "0"]
					output.write("\t".join(new_row) + "\n")
			# make sure we always start at 0 on a new contig and search for gaps
			print(f"detected gap at beginning of {contig}")
			last_ending = 0
			size_of_gap = start
			num_empty_entries = size_of_gap / interval_size
			print(f"Would create {num_empty_entries} 0 weight entries here")
			for s in range(last_ending, start, interval_size):
				print(f"\t{contig} {s} - {s + interval_size}")
				new_row = [contig, str(s), str(s + interval_size), ".", "0"]
				output.write("\t".join(new_row) + "\n")
			
		last_contig = contig
		last_ending = ending
		# now that we've inserted any fake rows necessary to fill gaps, write the current one
		output.write("\t".join(row) + "\n")
		
		
	#end for loop
	# check for anything leftover at the end of the last contig
	if last_contig in contig_sizes:
		# we'll want to round up the listed ending to the next block in case we
		# haven't reached it yet
		new_ending_block = int(math.ceil(contig_sizes[last_contig] / float(interval_size))) * interval_size
		print(f"Ending detected for contig: {last_contig} rounded up to block {new_ending_block}.  Last block written was {last_ending}")
		for s in range(last_ending, new_ending_block, interval_size):
			print(f"\t{last_contig} {s} - {s + interval_size}")
			new_row = [last_contig, str(s), str(s + interval_size), ".", "0"]
			output.write("\t".join(new_row) + "\n")
	print("Processing done")


