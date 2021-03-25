import sys
import gzip
import itertools

# Add new header for MAX_AS_VQSLOD

with gzip.open(sys.argv[1], 'rt') as file1:
    for line in file1:
        line = line.strip()

        if "##" in line:
            print(line)
            continue
        
        if "#CHROM" in line:
            print('##INFO=<ID=MAX_AS_VQSLOD,Number=1,Type=Float,Description="Maximum of AS_VQSLOD scores">')
            print(line)
            continue

        parts = line.split("\t")
        
        if (parts[6] == "ExcessHet"):
            continue
            
        info = parts[7]    
        d = dict([ tuple(x.split("=")) for x in info.split(";") if "=" in x]) 
        
        if "," in d['AS_VQSLOD']:
            pieces = [x for x in d['AS_VQSLOD'].split(",") if (x != "." and x != "NaN") ]
            
            if (len(pieces) == 1):
                m = pieces[0]
            elif (len(pieces) == 0):
                m = "."
            else:
                m = max([float(x) for x in pieces])
        else:
            m = d['AS_VQSLOD']

        parts[7] = f"{info};MAX_AS_VQSLOD={m}"
        print("\t".join(parts))
