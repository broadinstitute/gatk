import argparse
import sys

# A constant for comparing two floating point numbers
EPSILON = 0.00001

def compare_gvs_vcf_with_vds_vcf(gvs_vcf, vds_vcf, skip_gvs_filtered_lines):
    with open(gvs_vcf) as gvs:
        with open(vds_vcf) as vds:
            gvs_line = skip_header(gvs)
            vds_line = skip_header(vds)

            # Parse out the sample orders
            gvs_sample_to_column = get_sample_to_column(gvs_line)
            vds_sample_to_column = get_sample_to_column(vds_line)
            if (len(gvs_sample_to_column) != len(vds_sample_to_column)):
                print(f"DIFF: VCFs contain different number of samples: {len(gvs_sample_to_column)} vs {len(vds_sample_to_column)}")

            gvs_samples = sorted(gvs_sample_to_column.keys())
            vds_samples = sorted(vds_sample_to_column.keys())

            if (gvs_samples != vds_samples):
                print(f"DIFF: samples differ between VCFs:\n{gvs_line} vs \n{vds_line}")

            gvs_line = gvs.readline().rstrip()
            vds_line = vds.readline().rstrip()

            # sync up lines by position between the two VCFs
            count = 0
            while True:
                count += 1
                if skip_gvs_filtered_lines:
                    gvs_line = skip_filtered_lines(gvs_line, gvs)

                gvs_tokens = gvs_line.split()
                vds_tokens = vds_line.split()

                if (len(gvs_tokens) != len(vds_tokens)):
                    print(f"DIFF: different number of tokens in VCF line:\n{gvs_line} vs \n{vds_line}")
                    sys.exit(1)

                if (gvs_tokens[0] != vds_tokens[0]):
                    # TODO - need to deal with rolling over chromosomes between files
                    print(f"DIFF: CHROM differs between VCF line:\n{gvs_tokens[0]} vs {vds_tokens[0]}")
                    sys.exit(1)

                if (gvs_tokens[1] != vds_tokens[1]):
                    print(f"DIFF: POS differs between VCF line:\n{gvs_tokens[1]} vs {vds_tokens[1]}")
                    sys.exit(1)

                locus = gvs_tokens[0] + ":" + gvs_tokens[1]
                # print(f"{locus}")
                # if (int(gvs_tokens[1]) < 280726):
                #     gvs_line = gvs.readline().rstrip()
                #     vds_line = vds.readline().rstrip()
                #     continue

                if not args.verbose:
                    print('.', end='', flush=True)

                col_diff(locus, "ID", gvs_tokens[2], vds_tokens[2])
                col_diff(locus, "REF", gvs_tokens[3], vds_tokens[3])

                gvs_tokens[4], gvs_old_allele_index_to_new = reorder_gvs_alt_alleles(gvs_tokens[4])
                if (gvs_tokens[4] != vds_tokens[4]):
                    print(f"DIFF: ALT differs between VCF line:\n{gvs_tokens[4]} vs {vds_tokens[4]}")
                    if not args.report_all_diffs:
                        sys.exit(1)

                col_diff(locus, "QUAL", gvs_tokens[5], vds_tokens[5])

                # FILTER Field
                if vds_tokens[6] == "PASS":
                    vds_tokens[6] = "."
                col_diff(locus, "FILTER", gvs_tokens[6], vds_tokens[6])

                # INFO Field.
                # The info field is tricky. We only compare AC, AN, and AF. But VDS stores these differently
                # It includes the values for these metrics for the ref allele too.
                # Plus we have to deal with the ordering of alleles being different between GVS and VDS.
                if args.verbose:
                    print(f"\nLocus {locus}.  gvs INFO: {gvs_tokens[7]} vds INFO: {vds_tokens[7]}")

                gvs_info_fields = get_info_fields(gvs_tokens[7])
                vds_info_fields = get_info_fields(vds_tokens[7])

                if "AN" in gvs_info_fields:
                    info_diff(locus, gvs_info_fields, vds_info_fields, "AN", gvs_tokens[7], vds_tokens[7])
                    # Calculate the (expected) AC from the gvs VCF
                    gvs_an = int(gvs_info_fields["AN"])
                    reordered_gvs_ac = reorder_gvs_ac(gvs_info_fields["AC"], gvs_old_allele_index_to_new)
                    reordered_gvs_acs = reordered_gvs_ac.split(",")
                    alt_allele_ac = 0
                    for gvs_ac in reordered_gvs_acs:
                        alt_allele_ac += int(gvs_ac)
                    gvs_afs = []
                    gvs_afs.append(float(gvs_an - alt_allele_ac) / gvs_an)
                    new_gvs_ac = str(gvs_an - alt_allele_ac)        # AC for the ref
                    for gvs_ac in reordered_gvs_acs:
                        new_gvs_ac += "," + gvs_ac
                        gvs_afs.append(float(gvs_ac) / gvs_an)

                    vds_ac = vds_info_fields["AC"]
                    if (new_gvs_ac != vds_ac):
                        print(f"DIFF: Locus {locus}, INFO Field Value for AC differs between gvs and vds: {new_gvs_ac} vs {vds_ac}\nGVS INFO: {gvs_tokens[7]}\nVDS INFO: {vds_tokens[7]}")
                        if not args.report_all_diffs:
                            sys.exit(1)

                    vds_afs_string = vds_info_fields["AF"]
                    vds_afs = vds_afs_string.split(",")
                    if (len(gvs_afs) != len(vds_afs)):
                        print(f"GVS INFO AF {gvs_afs} has different number of elements than VDS INFO AF {vds_afs_string} alleles")
                        sys.exit(1)

                    for i in range(0, len(gvs_afs)):
                        vds_af = float(vds_afs[i])
                        if vds_af != 0.0:
                            diff = (gvs_afs[i] - vds_af) / vds_af

                        if abs(diff) > EPSILON:
                            print(f"DIFF: Locus {locus}, INFO Field Value for AF, element {i} differs between gvs and vds: {gvs_afs[i]} vs {vds_af}\nGVS INFO: {gvs_tokens[7]}\nVDS INFO: {vds_tokens[7]}")
                            if not args.report_all_diffs:
                                sys.exit(1)

                # FORMAT Field (it is just going to be different - we build a map of the fields to column for handling sample values)
                gvs_format_fields_to_column = get_format_fields_to_column(gvs_tokens[8])
                vds_format_fields_to_column = get_format_fields_to_column(vds_tokens[8])

                if args.verbose:
                    print(f"\nLocus {locus}.  gvs FORMAT: {gvs_tokens[8]} vds FORMAT: {vds_tokens[8]}")

                for sample in gvs_samples:
                    if args.verbose:
                        print(f"Sample {sample}.    gvs sample: {gvs_tokens[gvs_sample_to_column[sample]]}  vds_sample: {vds_tokens[vds_sample_to_column[sample]]}")

                    gvs_sample_genotypes = get_sample_genotype_fields(gvs_format_fields_to_column, gvs_tokens[gvs_sample_to_column[sample]])
                    vds_sample_genotypes = get_sample_genotype_fields(vds_format_fields_to_column, vds_tokens[vds_sample_to_column[sample]])

                    # NOTE: More complex rules for FT
                    # sample_key_diff(locus, sample, gvs_sample_genotypes, vds_sample_genotypes, "FT")

                    gvs_ft = "PASS"
                    if "FT" in gvs_sample_genotypes:
                        gvs_ft = gvs_sample_genotypes["FT"]
                        if gvs_ft == ".":
                            gvs_ft = "PASS"
                    vds_ft = vds_sample_genotypes["FT"]
                    if vds_ft == ".":
                        vds_ft = "PASS"
                    # For purposes of comparison a 'low_VQSLOD_INDEL' or 'low_VQSLOD_SNP' is recoded as 'FAIL' which is what the vds VCF uses.
                    if (gvs_ft == "low_VQSLOD_INDEL" or gvs_ft == "low_VQSLOD_SNP"):
                        gvs_ft = "FAIL"
                    if (gvs_ft != vds_ft):
                        print(f"\nDIFF: Locus {locus}, Sample {sample} Value for FT differs between gvs and vds: {gvs_ft} vs {vds_ft}")
                        print(f"Locus {locus}.  gvs FORMAT: {gvs_tokens[8]} vds FORMAT: {vds_tokens[8]}")
                        print(f"Sample {sample}.    gvs sample: {gvs_tokens[gvs_sample_to_column[sample]]}  vds_sample: {vds_tokens[vds_sample_to_column[sample]]}")
                        if not args.report_all_diffs:
                            sys.exit(1)

                    # Note on comparison here.
                    # Reorder GT string from GVS based on difference ordering of ALT alleles between GVS and VDS
                    # if VDS.LA and VDS.GT are both not "." (empty/missing)
                    # then generate NEW.VDS.GT using VDS.LGT and VDS.LA
                    # If the GT in the VDS generated VCF is "./."" AND VDS.FT == "FAIL" then VDS.GT = NEW.VDS.GT
                    # if VDS.GT != GVS.GT that's a difference (unless they are phased


                    # NOTE: More complex rules for GT
                    gvs_gt = reorder_gvs_gts(gvs_sample_genotypes["GT"], gvs_old_allele_index_to_new)
                    vds_la = vds_sample_genotypes["LA"]
                    vds_lgt = vds_sample_genotypes["LGT"]
                    vds_gt = vds_sample_genotypes["GT"]
                    new_vds_gt = vds_gt
                    if (vds_la != "." and vds_lgt != "."):
                        new_vds_gt = calculate_vds_gt_lgt(vds_lgt, vds_la)

                    if (vds_gt == "./." and vds_ft == "FAIL"):
                        # NOTE/HMM. In the VDS we no call a genotype if the FT is a FAIL, the original genotype is in LGT.
                        vds_gt = new_vds_gt
                    if not compare_gts(gvs_gt, vds_gt):
                    # if (gvs_gt != vds_gt):
                        print(f"\nDIFF: Locus {locus}, Sample {sample} Value for GT differs between gvs and vds: {gvs_gt} vs {vds_gt}")
                        print(f"Locus {locus}.  gvs FORMAT: {gvs_tokens[8]} vds FORMAT: {vds_tokens[8]}")
                        print(f"Sample {sample}.    gvs sample: {gvs_tokens[gvs_sample_to_column[sample]]}  vds_sample: {vds_tokens[vds_sample_to_column[sample]]}")
                        if not args.report_all_diffs:
                            sys.exit(1)

                    sample_key_diff(locus, sample, gvs_sample_genotypes, vds_sample_genotypes, "GQ")
                    sample_key_diff(locus, sample, gvs_sample_genotypes, vds_sample_genotypes, "RGQ")

                    ## AD
                    if "AD" in gvs_sample_genotypes:
                        gvs_ad = gvs_sample_genotypes["AD"]
                        vds_ad = vds_sample_genotypes["LAD"]
                        if (gvs_ad != vds_ad):
                            print(f"\nDIFF: Locus {locus}, Sample {sample} Value for differs between gvs key AD and vds key LAD: {gvs_ad} vs {vds_ad}")
                            print(f"Locus {locus}.  gvs FORMAT: {gvs_tokens[8]} vds FORMAT: {vds_tokens[8]}")
                            print(f"Sample {sample}.    gvs sample: {gvs_tokens[gvs_sample_to_column[sample]]}  vds_sample: {vds_tokens[vds_sample_to_column[sample]]}")
                            if not args.report_all_diffs:
                                sys.exit(1)

                gvs_line = gvs.readline().rstrip()
                vds_line = vds.readline().rstrip()

            exit(0)

def get_info_fields(info_string):
    """
    Method to split the fields from the INFO field into a dictionary
    :param info_string: string containing the INFO field
    :return: dictionary containing a map of key to value from the INFO field.
    """
    info_dict = {}
    tokens = info_string.split(";")
    for token in tokens:
        key_value = token.split("=")
        if len(key_value) != 2:
            print(f"INFO token {key_value} does not look right")
            sys.exit(1)
        info_dict[key_value[0]] = key_value[1]
    return info_dict

def compare_gts(gvs_gt, vds_gt):
    """
    Method to compare to GTs. Allows for equality between "0/1" and "1/0" for instance
    :param gvs_gt: GT field from the gvs VCF
    :param vds_gt: GT field from the vds VCF
    :return: true if the GT fields are the same
    """
    if (gvs_gt == vds_gt):
        return True

    tmp = gvs_gt.split("/")
    if len(tmp) != 2:
        print(f"GVS GT {gt} does not two elements - are these phased genotypes?? ('|' separator)")
        sys.exit(1)
    flipped_gvs_gt = tmp[1] + "/" + tmp[0]
    return flipped_gvs_gt == vds_gt

def calculate_vds_gt_lgt(lgt, la):
    """
    VDS stores the LGT field (which is the 'local genotypes')
    this represents the GT field as for ONLY the alleles available for the sample in question.
    We convert this to GT using LGT and the LA ('local alleles') field
    :param lgt: The LGT field from the vds VCF
    :param la: The LA field from the vds VCF
    :return: string containing the GT field
    """
    local_alleles = la.split(",")       # looks like "0,2" or "0,1,2"
    if len(local_alleles) != 2 and len(local_alleles) != 3:
        print(f"LA {la} does not have exactly two OR three elements!")
        sys.exit(1)
    local_gts = lgt.split("/")          # looks like "0/1"
    if len(local_gts) != 2:
        print(f"LGT {lgt} does not have exactly two elements!  - are these phased genotypes?? ('|' separator)")
        sys.exit(1)

    gt = local_alleles[int(local_gts[0])] + "/" + local_alleles[int(local_gts[1])]
    return gt

def reorder_gvs_ac(gvs_ac, gvs_old_allele_index_to_new):
    """
    Reorder the gvs's INFO AC field to be in the same order as the VDS's alleles
    :param gvs_ac: The AC field from the gvs VCF
    :param gvs_old_allele_index_to_new: A dictionary of old alt allele index to new
    :return: string containing the reordered AC for the gvs VCF (note that this now contains a value for REF too)
    """
    gvs_acs = gvs_ac.split(",")
    # Note that gvs_old_allele_index_to_new contains a mapping for the REF allele too, so one more expected than for ALT alleles
    if (len(gvs_acs) != len(gvs_old_allele_index_to_new.keys()) - 1):
        print(f"GVS INFO AC {gvs_ac} has different number of elements than GVS's alleles")
        sys.exit(1)
    reordered_gvs_acs = {}
    # Generate a new/better hash of old index to new (as ints)  TODO - probably should do it everywhere?
    old_index_to_new = {}
    for old in gvs_old_allele_index_to_new.keys():
        new_val = gvs_old_allele_index_to_new[old]
        old_val = int(old) - 1
        old_index_to_new[old_val] = int(new_val) - 1

    for i in range(0, len(gvs_acs)):
        reordered_gvs_acs[old_index_to_new[i]] = gvs_acs[i]

    reordered_gvs_acs_string = reordered_gvs_acs[0]
    for i in range(1, len(gvs_acs)):
        reordered_gvs_acs_string += "," + reordered_gvs_acs[i]

    return reordered_gvs_acs_string


def reorder_gvs_gts(gvs_gt, gvs_old_allele_index_to_new):
    """
    Since the ordering of the alt alleles in the vds VCF is different from that in the gvs VCF
    We need to reorder the genotype field correspondingly.
    This method does that:
    Remap genotypes for difference in ALT allele ordering between gvs VCF and vds VCF

    :param gvs_gt: The GT string from the gvs VCF (of the form "0/1")
    :param gvs_old_allele_index_to_new: A dictionary of the old allele index to the new
    (where old is that found in the gvs VCF and new is that in the vds VCF)
    :return: gvs_gt string representing the newly encoded GT (of the form "0/2")
    """
    if gvs_gt != "./.":
        gvs_gts = gvs_gt.split("/")
        if len(gvs_gts) != 2:
            print(f"GVS GT {gvs_gt} does not two elements - are these phased genotypes?? ('|' separator)")
            sys.exit(1)
        new_gvs_gts = []
        for i in range(0, len(gvs_gts)):
            gvs_gt = gvs_gts[i]
            if gvs_gt not in gvs_old_allele_index_to_new:
                print(f"Didn't find it!!")
                sys.exit(1)
            new_gvs_gts.append(gvs_old_allele_index_to_new[gvs_gt])
        gvs_gt = "/".join(new_gvs_gts)
    return gvs_gt

def reorder_gvs_alt_alleles(alt_allele_string):
    """
    The alt alleles in gvs are ordered differently than those in vds. gvs seems to be by first usage,
    vds are ordered alphabetically. Here we reorder (by simple sort) the alt alleles in the gvs record
    and generate a directory of old GT (1, 2, 3) to new (reordered) GT (2, 3, 1) for instance.

    :param alt_allele_string: The alt allele from the gvs VCF (of the form "A,ACTAA,ACT")
    :return: the reordered alt_allele_string (in the new order).
    :return: A map of the old alt allele position to the new position in the reordered alt_allele_string
    """
    gvs_alleles = alt_allele_string.split(",")
    if (len(gvs_alleles) == 0):
        print(f"No ALT ALLELES???")
        sys.exit(1)

    gvs_old_allele_index_to_new = {}
    gvs_old_allele_index_to_new["0"] = "0"      # Put in an entry for ref.
    if (len(gvs_alleles) > 1):
        old_alleles = gvs_alleles[:]        # copy the array
        gvs_alleles.sort()
        alt_allele_string = ",".join(gvs_alleles)

        # Make a dictionary of the OLD alt allele index to the new (SORTED)
        for i in range(0, len(old_alleles)):
            for j in range(0, len(gvs_alleles)):
                if old_alleles[i] == gvs_alleles[j]:
                    gvs_old_allele_index_to_new[str(i+1)] = str(j+1)
                    break
    else:
        gvs_old_allele_index_to_new["1"] = "1"
    return(alt_allele_string, gvs_old_allele_index_to_new)

def col_diff(locus, key, gvs_value, vds_value):
    """
    Simple method to compare two values and give context sensitive information if they differ
    :param locus: The locus (chr:pos) in the VCF
    :param key: The key (column) in the VCF fom which the fields are pulled
    :param gvs_value: The value from the gvs
    :param vds_value: The value from the vds
    :return: Nothing
    """
    if (gvs_value != vds_value):
        print(f"DIFF: Locus {locus}, {key} Value differs between gvs and vds: {gvs_value} vs {vds_value}")
        if not args.report_all_diffs:
            sys.exit(1)

def info_diff(locus, gvs_dict, vds_dict, key, gvs_info, vds_info):
    """
    Simple method to compare two values in the INFO field and give context sensitive information if they differ
    :param locus: The locus (chr:pos) in the VCF
    :param gvs_dict: A dictionary of the key:values in the gvs INFO field
    :param vds_dict: A dictionary of the key:values in the vds INFO field
    :param key: The key (column) in the INFO field dictionary from which the values will be pulled
    :param gvs_info: The INFO string from the gvs, used for logging
    :param vds_info: The INFO string from the vds, used for logging
    :return: Nothing
    """
    gvs_value = gvs_dict[key]
    vds_value = vds_dict[key]
    if (gvs_value != vds_value):
        print(f"DIFF: Locus {locus}, INFO Field Value for {key} differs between gvs and vds: {gvs_value} vs {vds_value}\nGVS INFO: {gvs_info}\nVDS INFO: {vds_info}")
        if not args.report_all_diffs:
            sys.exit(1)

def sample_key_diff(locus, sample, gvs_dict, vds_dict, key):
    """
    Simple method to compare two values in the SAMPLE field and give context sensitive information if they differ
    :param locus: The locus (chr:pos) in the VCF
    :param sample: The sample name, used for logging
    :param gvs_dict: A dictionary of the key:values in the gvs SAMPLE field
    :param vds_dict: A dictionary of the key:values in the vds SAMPLE field
    :param key: The key (column) in the SAMPLE field dictionary from which the values will be pulled
    :return: Nothing
    """
    gvs_value = gvs_dict[key]
    vds_value = vds_dict[key]
    if (gvs_value != vds_value):
        print(f"DIFF: Locus {locus}, Sample {sample} Value for {key} differs between gvs and vds:\n{gvs_value} vs {vds_value}")
        if not args.report_all_diffs:
            sys.exit(1)

def get_sample_genotype_fields(format_fields_to_column, sample_genotype_string):
    """
    A method to parse the fields out of the sample_genotype_string and return a map of them.
    :param format_fields_to_column: A map of the FORMAT field names to their column value.
    :param sample_genotype_string: The sample genotype string
    :return: a dictionary of format field to sample value (eg. GQ -> 0.2)
    """
    tokens = sample_genotype_string.split(":")
    sample_genotype_fields = {}
    for key, index in format_fields_to_column.items():
        value = "."         # Missing field
        if (index < len(tokens)):
            value = tokens[index]
        sample_genotype_fields[key] = value
    return sample_genotype_fields

def get_format_fields_to_column(token):
    """
    A method to parse the FORMAT field and return it a map of format field name to column index
    :param token: The FORMAT field from the VCF
    :return: a dictionary of the format field names to their index in the string
    """
    tokens = token.split(":")
    fields_to_column = {}
    for i in range(0, len(tokens)):
        fields_to_column[tokens[i]] = i
    return fields_to_column

def get_sample_to_column(fline):
    """
    A method to parse all of the sample columns from the VCF and return a map of sample to column index
    :param fline: The VCF header string
    :return: a dictionary of sample name to column index
    """
    tokens = fline.split()
    sample_to_column = {}
    for i in range(9, len(tokens)):
        sample_to_column[tokens[i]] = i
    return sample_to_column

def skip_filtered_lines(fline, fp):
    """
    A method to skip over any 'filtered' lines in the VCF. If the line has a filter applied, it will skip that line
    :param fline: The latest line (as a string) from the VCF
    :param fp: The file's file pointer.
    :return: The nextmost line from the file that has not been filtered.
    """
    ftokens = fline.split()
    filter = ftokens[6]
    while True:
        if filter == "." or filter == "PASS":
            return fline
        fline = fp.readline().rstrip()
        ftokens = fline.split()
        filter = ftokens[6]

def skip_header(fp):
    """
    A method to skip over the header in the VCF
    :param fp: The file's file pointer
    :return: The nextmost line that is not a header line
    """
    while (line := fp.readline().rstrip()):
        if not line.startswith("##"):
            return line


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='A tool to compare a GVS-generated VCF with one generated from VDS')

    parser.add_argument('--gvs_vcf', type=str, help='GVS-generated VCF', required=True)
    parser.add_argument('--vds_vcf', type=str, help='VDS-generated VCF', required=True)
    parser.add_argument('--skip_gvs_filtered_lines', action='store_true', help='If set, skip lines that are filtered in the GVS-generated VCF')
    parser.add_argument("--verbose", action="store_true", help="increase output verbosity")
    parser.add_argument("--report_all_diffs", action="store_true", help="Do NOT exit on the first difference found")
    args = parser.parse_args()

    compare_gvs_vcf_with_vds_vcf(args.gvs_vcf, args.vds_vcf, args.skip_gvs_filtered_lines)
