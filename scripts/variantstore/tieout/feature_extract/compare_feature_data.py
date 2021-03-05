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
            alts = parts[4]
            if ("*" in alts or loc in exclude_set):
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
    
    errors = [ key for key in ANNOTATIONS if not equals(i1, i2, key) ]
    if len(errors) > 0:
        error_fields = ",".join(errors)
        print(f"DIFF on {error_fields} {e1['chrom']}:{e1['pos']} for {e1['ref']}/{e1['alt']} and {e2['ref']}/{e2['alt']}")
        print(i1)
        print(i2)
        print("--------------")
                
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
        for key in ['chrom','pos','id', 'ref', 'alt']:
            if not equals(e1, e2, key):
                log_difference(key, e1, e2) 

        compare_features(e1, e2)
        
        # compare the minimized version of ref/alt
#        compare_alts(e1['alt'], e2['alt'])

#        print(f"Success comparing {e1['pos']}")
#        sys.exit(0)
        
        
        lines = lines + 1
        if (lines % 100000 == 0):
            print(f"Compared {lines} positions")
        
        
print(f"Compared {lines} positions")