import argparse
import pandas as pd

location_offset = 1000000000000

argParser = argparse.ArgumentParser()

argParser.add_argument("-s", "--bin-size", help="size of the bins in kb", type=int, default=1)
argParser.add_argument("-i", "--infile", help="bed input file", type=str, default="49K_vet_weight_1k.csv")
argParser.add_argument("-o", "--outfile", help="bed output file", type=str, default="gvs_vet_weights_1kb.bed")


args = argParser.parse_args()
print("args=%s" % args)

binsize_kb = args.bin_size
infile = args.infile
outfile = args.outfile

print(f"Input file is {infile}")

w = pd.read_csv(infile, dtype={'bin': int, 'entries':int})

w['contig'] = "chr" + (w['bin'].astype(int) / location_offset).astype(int).astype(str).str.replace("23","X").replace("24","Y")
w['start_position'] = w['bin'].astype(int) - (w['bin'] / location_offset).astype(int) * location_offset
w['end_position'] = w['start_position'] + binsize_kb*1000
w['name'] = "."
o = w[['contig', 'start_position','end_position', 'name', 'entries' ]]

o.to_csv(outfile,sep='\t',index=False,header=False)