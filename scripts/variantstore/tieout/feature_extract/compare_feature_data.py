import sys
import gzip
import itertools

ANNOTATIONS = ["AS_FS","AS_MQ","AS_MQRankSum","AS_QD","AS_ReadPosRankSum","AS_SOR"]

def get_next_line(i):
    for line in i:
        if (line.startswith("##")):
            pass;
        else:
            parts = line.strip().split("\t")
            loc = f"{parts[0]}:{parts[1]}"
            if (loc in exclude_set):
#                print(f"Skipping {loc}")
                pass;
            else:
                return line;

def parseline(e, header):
    data = {}
    
    parts = e.strip().split("\t")

    data['chrom'] = parts[0]
    data['pos'] = int(parts[1])
    data['id'] = parts[2]
    ref= parts[3]
    data['orig_alt'] = parts[4]
    
    alts = [x for x in parts[4].split(",")]
    
    ## and now minimize the ref and alt (with * allele being excempted)
    done = False
    while (not done and len(ref) != 1):
        if (all(ref[-1] == alt[-1] for alt in alts if alt != '*')):
            ref = ref[:-1]
            alts = [alt[:-1] if alt != '*' else '*' for alt in alts]
        else:
            done = True
    
    data['ref'] = ref
    data['alt'] = ",".join(alts)
    
    data['filter'] = parts[6]
    
    d = dict([ tuple(x.split("=")) for x in parts[7].split(";") ]) 
    data['info_data'] = {key: d[key] for key in d if key in ANNOTATIONS}
    
    return data;

def equals(e1, e2, key):
    return (key in e1 and key in e2 and e1[key] == e2[key])

def compare_float(e1, e2, key, tolerance, location):
    # compare directly first, also handles '.' case
    # TODO: deal with multi allelics properly
    s1 = e1[key].split(",")[0] if "," in e1[key] else e1[key]
    s2 = e2[key]

    if (s1 != s2):
        if ("." == s1 or "." == s2):
            print(f"DIFF on {key} at {location} with values of {s1} and {s2}")
            print(f"{e1}")
            print(f"{e2}")
            print("--------------")

        else:
            v1 = float(s1)
            v2 = float(s2)

            delta = abs(v2 - v1)
            if delta > tolerance:
                print(f"DIFF on {key} at {location} of {delta}")
                print(f"{e1}")
                print(f"{e2}")
                print("--------------")
                

# need to sort, order doesn't matter (at this point)
def compare_alts(e1, e2):
    p1 = [x for x in e1.split(",") if x != '*']
    p1.sort()
    
    p2 = [x for x in e2.split(",") if x != '*']
    p2.sort()
    
    s1 = ",".join(p1)
    s2 = ",".join(p2)

    if (s1 != s2):
        print(f"DIFF on ALTS")
        print(f"{s1}")
        print(f"{s2}")

def get_gt_indexes(gt):
    delim = "|" if "|" in gt else "/"
    
    gt1 = gt.split(delim)[0]
    gt2 = gt.split(delim)[1]
    return (gt1, gt2)

def get_gt_alleles(gt, ref, alts):
    alleles = [ref] + alts.split(",")
    
    (gt1, gt2) = get_gt_indexes(gt)
        
    a1 = alleles[int(gt1)] if gt1 != "." else "."        
    a2 = alleles[int(gt2)] if gt2 != "." else "."

    # TODO: for now, ignore phasing...
    return [a1,a2]

def log_difference(key, e1, e2, sample_id = None):
    if sample_id:
        sd1 = e1['sample_data']
        sd2 = e2['sample_data']
        a1 = get_gt_alleles(sd1[sample_id]['GT'], e1['ref'], e1['alt'])
        a2 = get_gt_alleles(sd2[sample_id]['GT'], e2['ref'], e2['alt'])
                
        print(f"DIFF on {key} for {sample_id} at {e1['chrom']}:{e1['pos']} with {e1['alt']} and {e2['alt']}")
        print(f"{a1} vs {a2}")
        print(sd1[sample_id])
        print(sd2[sample_id])
    else:
        print(f"DIFF on {key} {e1['chrom']}:{e1['pos']} with {e1['alt']} and {e2['alt']}")
        print(e1)
        print(e2)

    print("--------------")
        
def compare_features(e1, e2):
    i1 = e1['info_data']
    i2 = e2['info_data']
    
    # TODO: deal with multi-allelics eventually...
    for key in ANNOTATIONS:
        compare_float(i1, i2, key, 0.0, f"{e1['chrom']}:{e1['pos']}")
        
def unroll_interval_range(r):
    (chrom, range_string) = r.split(":")
    (start, end) = range_string.split("-")
    return [ f"{chrom}:{x}" for x in range(int(start), int(end)+1) ]
    
vcf_file_1 = sys.argv[1]
vcf_file_2 = sys.argv[2]

exclude_list = []
if (len(sys.argv) == 4):
    with open(sys.argv[3]) as f:
        for x in f.readlines():
            if "-" not in x:
                exclude_list.append(x)
            else:
                exclude_list.extend(unroll_interval_range(x.strip()))
exclude_set = set(exclude_list)


print(f"Excluding {len(exclude_set)} loci")

lines = 0
with gzip.open(vcf_file_1, 'rt') as file1, gzip.open(vcf_file_2, 'rt') as file2:

    while True:
        line1 = get_next_line(file1)
        line2 = get_next_line(file2)

        if (line1 == None and line2 == None):
            break

        # we don't care about samples, so skip the sample-header line as well
        if (line1.startswith("#")):
            continue;
            
        # parse out data
        e1 = parseline(line1, None)
        e2 = parseline(line2, None)

        # skipping multi-alleleics for now
        if (len(e1['alt'].replace(",","").replace("*","")) > 1):
            print(f"Skipping {e1['chrom']}:{e1['pos']} since it is multiallelic")
            continue;

        while (e1['pos'] != e2['pos']):
            if (e1['pos'] < e2['pos']):
                print(f"DIFF on position {e1['chrom']}:{e1['pos']} is MISSING in file 2.  Advancing")
                print("--------------")
                line1 = get_next_line(file1)    
                e1 = parseline(line1, None)
            else:
                print(f"DIFF on position {e2['chrom']}:{e2['pos']} is MISSING in file 1.  Advancing")
                print("--------------")
                line2 = get_next_line(file2)
                e2 = parseline(line2, None)

        
        # do the comparison of exact matches at the position level
        for key in ['chrom','pos','id', 'ref']:
            if not equals(e1, e2, key):
                log_difference(key, e1, e2) 

        compare_features(e1, e2)
        
        # compare the minimized version of ref/alt
        compare_alts(e1['alt'], e2['alt'])

        print(f"Success comparing {e1['pos']}")
#        sys.exit(0)
        
        
        lines = lines + 1
        if (lines % 100000 == 0):
            print(f"Compared {lines} positions")
        
        
print(f"Compared {lines} positions")