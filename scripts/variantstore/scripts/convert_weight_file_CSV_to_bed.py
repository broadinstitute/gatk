import argparse
import pandas as pd
import csv

location_offset = 1000000000000
hg38_contig_map = {1: 'chr1', 2: 'chr2', 3: 'chr3', 4: 'chr4', 5: 'chr5', 6: 'chr6', 7: 'chr7', 8: 'chr8', 9: 'chr9', 10: 'chr10', 11: 'chr11', 12: 'chr12', 13: 'chr13', 14: 'chr14', 15: 'chr15', 16: 'chr16', 17: 'chr17', 18: 'chr18', 19: 'chr19', 20: 'chr20', 21: 'chr21', 22: 'chr22', 23: 'chrX', 24: 'chrY'}

argParser = argparse.ArgumentParser()

argParser.add_argument("-s", "--bin-size", help="size of the bins in kb", type=int, default=1)
argParser.add_argument("-i", "--infile", help="bed input file", type=str, default="49K_vet_weight_1k.csv")
argParser.add_argument("-o", "--outfile", help="bed output file", type=str, default="gvs_vet_weights_1kb.bed")
argParser.add_argument("-m", "--mapping-file", help="chromosome to integer mapping file", type=str)


args = argParser.parse_args()
print("args=%s" % args)

binsize_kb = args.bin_size
infile = args.infile
outfile = args.outfile

# Function to parse the contig mappings file
def load_contig_mapping(mapping_file):
    contig_map = {}
    with open(mapping_file, mode='r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) == 2:
                contig_map[int(row[1])] = row[0]
            else:
                raise ValueError(f"Invalid line format: {row}")
    return contig_map

# Load the contig mapping file if provided, else use standard hg38 mappings
if args.mapping_file:
    contig_map = load_contig_mapping(args.mappingfile)
else:
    contig_map = hg38_contig_map

w = pd.read_csv(infile, dtype={'bin': int, 'entries':int})

w['contig'] = (w['bin'].astype(int) / location_offset).astype(int).map(contig_map)
w['start_position'] = w['bin'].astype(int) - (w['bin'] / location_offset).astype(int) * location_offset
w['end_position'] = w['start_position'] + binsize_kb*1000
w['name'] = "."
o = w[['contig', 'start_position','end_position', 'name', 'entries' ]]

o.to_csv(outfile,sep='\t',index=False,header=False)