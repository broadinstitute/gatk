import sys
import gzip

# Add new header for MAX_AS_VQS_SENS and MAX_AS_VQSLOD

with gzip.open(sys.argv[1], 'rt') as file1:
    for line in file1:
        line = line.strip()

        if "##" in line:
            print(line)
            continue
        
        if "#CHROM" in line:
            print('##INFO=<ID=MAX_AS_VQS_SENS,Number=1,Type=Float,Description="Maximum of AS_VQS_SENS scores">')
            print('##INFO=<ID=MAX_AS_VQSLOD,Number=1,Type=Float,Description="Maximum of AS_VQSLOD scores">')
            print(line)
            continue

        parts = line.split("\t")
        
        # strip out hard filtered sites, so vcfeval can use "all-records" to plot the ROC curves
        if ("ExcessHet" in parts[6] or "LowQual" in parts[6] or "NO_HQ_GENOTYPES" in parts[6]):
            continue

        info = parts[7]    
        d = dict([ tuple(x.split("=")) for x in info.split(";") if "=" in x]) 

        format_key = [x for x in parts[8].split(":")]
        sample_data = dict(zip(format_key, parts[9].split(":")))

        gt = sample_data['GT']
        
        if gt == "0/0" or gt == "./.":
            continue
            
        if 'FT' in sample_data:
            ft = sample_data['FT']

            # if there is a non-passing FT value
            if not (ft == "PASS" or ft == "."):
                
                # overwrite FILTER if it was PASS or "."
                if (parts[6] == "PASS" or parts[6] == "."):
                    parts[6] = ft
                # otherwise append it to the end
                else:
                    parts[6] = parts[6] + "," + ft

        if "AS_VQS_SENS" in d:
            if "," in d['AS_VQS_SENS']:
                pieces = [x for x in d['AS_VQS_SENS'].split(",") if (x != "." and x != "NaN") ]

                if (len(pieces) == 1):
                    m = pieces[0]
                elif (len(pieces) == 0):
                    m = "."
                else:
                    m = max([float(x) for x in pieces])
            else:
                m = d['AS_VQS_SENS']
            parts[7] = f"{info};MAX_AS_VQS_SENS={m}"
        elif "AS_VQSLOD" in d:
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
        else:
            sys.exit(f"Can find neither AS_VQS_SENS nor AS_VQSLOD in {line}")

        print("\t".join(parts))
