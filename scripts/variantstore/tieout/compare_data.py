import sys

def get_next_line(i):
    for line in i:
        if (line.startswith("##")):
            pass;
        else:
            parts = line.strip().split("\t")
            loc = f"{parts[0]}:{parts[1]}"
            if (loc in exclude_list):
#                print(f"Skipping {loc}")
                pass;
            else:
                return line;

def compare_headers(header1, header2):
    a = set(header1.split("\t"))
    b = set(header2.split("\t"))
    diff = a.symmetric_difference(b)
    
    if (len(diff) == 0):
        print("Headers match, including all samples...")
    else:
        print(f"DIFF: headers are different! {a} vs {b}")
        sys.exit(1)

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

    format_key = [x for x in parts[8].split(":")]

    samples = header.split("\t")[9:]
    sample_data = [ dict(zip(format_key, x.split(":"))) for x in parts[9:] ]
    
    data['sample_data'] =  dict(zip( samples, sample_data))

    return data;

def equals(e1, e2, key):
    return (key in e1 and key in e2 and e1[key] == e2[key])

def equals_int(e1, e2, key, tolerance):
    if (key in e1 and key in e2):
        v1 = int(e1[key]) if e1[key] != "." else -1
        v2 = int(e2[key]) if e2[key] != "." else -1
        
        return (abs(v2-v1) <= tolerance)
    else:
        return False

def compare_float(e1, e2, key, tolerance):
    # compare directly first, also handles '.' case
    s1 = e1[key]
    s2 = e2[key]

    if (s1 != s2):
        if ("." in s1 or "." in s2):
            print(f"DIFF on {key} with values of {e1} and {e2}")

        else:
            v1 = float(s1)
            v2 = float(s2)

            delta = abs(v2 - v1)
            if delta > tolerance:
                print(f"DIFF on {key} of {delta}")
                print(f"{e1}")
                print(f"{e2}")

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
        
def get_pl_for_gt(gt,pl):
    (gt1, gt2) = get_gt_indexes(gt)

    # No PLs for no-calls
    if (gt1 == "." or gt2 == "."):
        return None

    # from the VCF spec:  for P=2, the index of the genotype “a/b”, where a ≤ b, is b(b + 1)/2 + a
    a = min(int(gt1), int(gt2))
    b = max(int(gt1), int(gt2))
    i = int(b * (b+1)/2 + a)
    
#    print(f"Calculated PL index {i} for {gt} with a: {a} and b: {b}")

    return int(pl.split(",")[i])
        
def compare_sample_data(e1, e2):
    sd1 = e1['sample_data']
    sd2 = e2['sample_data']
    
    if (len(sd1) != len(sd2)):
        print(f"DIFF on length of sample data {len(sd1)} and {len(sd2)}")
        print(f"{e1}")
        print(f"{e2}")
 
    for sample_id in sd1.keys():
        # if either has a GQ... compare it!
        if 'GQ' in sd1[sample_id] or 'GQ' in sd2[sample_id]:
            if not equals_int(sd1[sample_id], sd2[sample_id], 'GQ', 2):
                log_difference('GQ', e1, e2, sample_id) 

        # if both have RGQ compare it (unlikely)
        if 'RGQ' in sd1[sample_id] and 'RGQ' in sd2[sample_id]:
            if not equals(sd1[sample_id], sd2[sample_id], 'RGQ'):
                log_difference('RGQ', e1, e2, sample_id) 

        # if first (older) has PLs, and seccond (bq) has RGQ... compare!
        elif 'PL' in sd1[sample_id] and 'RGQ' in sd2[sample_id]:
            rgq1 = sd1[sample_id]['PL'].split(",")[0]
            rgq2 = sd2[sample_id]['RGQ']
            
            if (rgq1 != rgq2):
                log_difference('RGQ', e1, e2, sample_id) 
            
        # compare genotypes based on actual alleles (since order of alts might differ)
        a1 = get_gt_alleles(sd1[sample_id]['GT'], e1['ref'], e1['alt'])
        a2 = get_gt_alleles(sd2[sample_id]['GT'], e2['ref'], e2['alt'])
        
        if (set(a1) != set(a2)):
            
            # TODO: hack to work around myriad of errors where overlapping reference blocks incorrectly report ./. in classic pipeline
            if ( sd1[sample_id]['GT'] == "./." and sd2[sample_id]['GT'] == "0/0"):
                return

            # If the genotypes are different BUT they have the same PL, they are effectively eqivalent.  

            if 'PL' in sd1[sample_id] and 'PL' in sd2[sample_id]:
                pl1 = get_pl_for_gt(sd1[sample_id]['GT'], sd1[sample_id]['PL'])
                pl2 = get_pl_for_gt(sd2[sample_id]['GT'], sd2[sample_id]['PL'])

                # print(f"comparing pl1 vs pl2: {pl1} vs {pl2}")
                if ( pl1 == pl2 ):
                    return

            # special case where WARP drops PLs, we accept both being GQ0 as equivalent
            if 'PL' not in sd1[sample_id] and int(sd1[sample_id]['GQ']) == 0 and int(sd2[sample_id]['GQ']) == 0:
                return

            log_difference('Genotypes', e1, e2, sample_id) 

    
vcf_file_1 = sys.argv[1]
vcf_file_2 = sys.argv[2]

exclude_list = []
if (len(sys.argv) == 4):
    with open(sys.argv[3]) as f:
        exclude_list = [x.strip() for x in f.readlines()]

lines = 0
with open(vcf_file_1) as file1, open(vcf_file_2) as file2:

    while True:
        line1 = get_next_line(file1)
        line2 = get_next_line(file2)

        if (line1 == None and line2 == None):
            break

        # get headers and compare (may have different sample order, that's ok)
        if (line1.startswith("#")):
            header1 = line1.strip()
            header2 = line2.strip()
            compare_headers(header1, header2)
            continue;
            
        # parse out data
        e1 = parseline(line1, header1)
        e2 = parseline(line2, header2)

        while (e1['pos'] != e2['pos']):
            if (e1['pos'] < e2['pos']):
                print(f"DIFF on position {e1['chrom']}:{e1['pos']} is MISSING in file 2.  Advancing")
                print("--------------")
                line1 = get_next_line(file1)    
                e1 = parseline(line1, header1)
            else:
                print(f"DIFF on position {e2['chrom']}:{e2['pos']} is MISSING in file 1.  Advancing")
                print("--------------")
                line2 = get_next_line(file2)
                e2 = parseline(line2, header2)

        # do the comparison of exact matches at the position level
        for key in ['chrom','pos','id', 'ref']:
            if not equals(e1, e2, key):
                log_difference(key, e1, e2) 

        # TODO: temporary until we decide what to do with spanning deletions
        if ('*' in e1['alt']):
#            print(f"Dropping {e1['chrom']}:{e1['pos']} due to * allele")
            continue
            
        # compare the minimized version of ref/alt
        compare_alts(e1['alt'], e2['alt'])

        # compare sample level data
        compare_sample_data(e1, e2)
        
        lines = lines + 1
        
        
print(f"Compared {lines} positions")
