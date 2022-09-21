import argparse
import sys
import numpy

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

                if (gvs_tokens[2] != vds_tokens[2]):
                    print(f"DIFF: ID differs between VCF line:\n{gvs_tokens[2]} vs {vds_tokens[2]}")
                    if not args.report_all_diffs:
                        sys.exit(1)

                if (gvs_tokens[3] != vds_tokens[3]):
                    print(f"DIFF: REF differs between VCF line:\n{gvs_tokens[3]} vs {vds_tokens[3]}")
                    if not args.report_all_diffs:
                        sys.exit(1)

                gvs_tokens[4], gvs_old_allele_index_to_new = reorder_gvs_alt_alleles(gvs_tokens[4])
                # print(f"-> {gvs_old_allele_index_to_new}")
                if (gvs_tokens[4] != vds_tokens[4]):
                    print(f"DIFF: ALT differs between VCF line:\n{gvs_tokens[4]} vs {vds_tokens[4]}")
                    if not args.report_all_diffs:
                        sys.exit(1)

                if (gvs_tokens[5] != vds_tokens[5]):
                    print(f"DIFF: QUAL differs between VCF line:\n{gvs_tokens[5]} vs {vds_tokens[5]}")
                    if not args.report_all_diffs:
                        sys.exit(1)

                if gvs_tokens[6] == "PASS":
                    gvs_tokens[6] = "."
                if vds_tokens[6] == "PASS":
                    vds_tokens[6] = "."

                if (gvs_tokens[6] != vds_tokens[6]):
                    print(f"DIFF: FILTER differs between VCF line:\n{gvs_tokens[6]} vs {vds_tokens[6]}")
                    if not args.report_all_diffs:
                        sys.exit(1)

                # Deal with differences in the INFO field
                if (gvs_tokens[7] != vds_tokens[7]):
                    print(f"\nDIFF: Locus {locus}, INFO differs between gvs and vds:\n{gvs_tokens[7]}\n vs \n{vds_tokens[7]}")
                    # if not args.report_all_diffs:
                    #     sys.exit(1)

                gvs_info_fields = get_info_fields(gvs_tokens[7])
                vds_info_fields = get_info_fields(vds_tokens[7])
                # TODO - handle new fields?

                if "AN" in gvs_info_fields:
                    info_diff(gvs_info_fields, vds_info_fields, "AN", locus, gvs_tokens[7], vds_tokens[7])
                    # Calculate the (expected) AC from the gvs VCF
                    gvs_an = int(gvs_info_fields["AN"])
                    reordered_gvs_ac = reorder_gvs_ac(gvs_info_fields["AC"], gvs_old_allele_index_to_new)
                    reordered_gvs_acs = reordered_gvs_ac.split(",")
                    alt_allele_ac = 0
                    for gvs_ac in reordered_gvs_acs:
                        alt_allele_ac += int(gvs_ac)
                    new_gvs_ac = str(gvs_an - alt_allele_ac)        # AC for the ref
                    for gvs_ac in reordered_gvs_acs:
                        new_gvs_ac += "," + gvs_ac

                    # new_gvs_af = gvs_ac / gvs_an
##                    new_gvs_ac = str(gvs_an - gvs_ac) + "," + str(gvs_ac)
                    # ref_allele_af = '{:0.6f}'.format((float) (gvs_an - gvs_ac) / gvs_an)
                    # new_vcs_ac = str(((int) vcs_an) - ((int) vcs_ac)) + "," + vcs_ac
                    # print(f"{ref_allele_af}")

                    vds_ac = vds_info_fields["AC"]
                    if (new_gvs_ac != vds_ac):
                        print(f"DIFF: Locus {locus}, INFO Field Value for AC differs between gvs and vds: {new_gvs_ac} vs {vds_ac}\nGVS INFO: {gvs_tokens[7]}\nVDS INFO: {vds_tokens[7]}")
                        if not args.report_all_diffs:
                            sys.exit(1)

                    # vds_af = vds_info_fields["AF"]
                    # if (new_gvs_af != vds_af):
                    #     print(f"DIFF: Locus {locus}, INFO Field Value for AF differs between gvs and vds: {new_gvs_af} vs {vds_af}\nGVS INFO: {gvs_tokens[7]}\nVDS INFO: {vds_tokens[7]}")
                    #     if not args.report_all_diffs:
                    #         sys.exit(1)

                    # info_diff(gvs_info_fields, vds_info_fields, "AF", locus, gvs_tokens[7], vds_tokens[7])

                # FORMAT is just going to be different
                gvs_format_fields_to_column = get_format_fields_to_column(gvs_tokens[8])
                vds_format_fields_to_column = get_format_fields_to_column(vds_tokens[8])

                if args.verbose:
                    print(f"\nLocus {locus}.  gvs FORMAT: {gvs_tokens[8]} vds FORMAT: {vds_tokens[8]}")
                else:
                    print('.', end='')

                for sample in gvs_samples:
                    if args.verbose:
                        print(f"Sample {sample}.    gvs sample: {gvs_tokens[gvs_sample_to_column[sample]]}  vds_sample: {vds_tokens[vds_sample_to_column[sample]]}")

                    gvs_sample_genotypes = get_sample_genotype_fields(gvs_format_fields_to_column, gvs_tokens[gvs_sample_to_column[sample]])
                    vds_sample_genotypes = get_sample_genotype_fields(vds_format_fields_to_column, vds_tokens[vds_sample_to_column[sample]])

                    # NOTE: More complex rules for FT
                    # simple_diff(gvs_sample_genotypes, vds_sample_genotypes, "FT", locus, sample)
                    gvs_ft = "PASS"
                    if "FT" in gvs_sample_genotypes:
                        # TODO - NOTE there are gvs records without a FT field
                        gvs_ft = gvs_sample_genotypes["FT"]
                    vds_ft = vds_sample_genotypes["FT"]
                    # TODO - this should be undone once the VDS-generated VCF is updated to have an explicity 'PASS'
                    if vds_ft == ".":
                        vds_ft = "PASS"
                    if gvs_ft == ".":
                        gvs_ft = "PASS"
                    # END TODO
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

                    simple_diff(gvs_sample_genotypes, vds_sample_genotypes, "GQ", locus, sample)
                    simple_diff(gvs_sample_genotypes, vds_sample_genotypes, "RGQ", locus, sample)

                    ## AD
                    if "AD" in gvs_sample_genotypes:        # TODO - did this go away???
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
                # if count > 10000:
                # sys.exit(1)

            exit(0)

def get_info_fields(info_string):
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
    if (gvs_gt == vds_gt):
        return True

    tmp = gvs_gt.split("/")
    if len(tmp) != 2:
        print(f"GVS GT {gt} does not two elements - are these phased genotypes?? ('|' separator)")
        sys.exit(1)
    flipped_gvs_gt = tmp[1] + "/" + tmp[0]
    return flipped_gvs_gt == vds_gt

def calculate_vds_gt_lgt(lgt, la):
    # VDS stores the LGT field (which is the 'local genotypes') this represents the GT field as for ONLY the alleles available for the sample in question.
    # We convert this to GT using LGT and the LA ('local alleles') field
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
    # Reorder the gvs's INFO AC field to be in the same order as the VDS's alleles
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
    # Remap genotypes for difference in ALT allele ordering between gvs VCF and vds VCF
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
    # So, the alt alleles in gvs are ordered differently than those in vds. gvs seems to be by first usage,
    # vds are ordered alphabetically. Here we reorder (by simple sort) the alt alleles in the gvs record
    # and generate a directory of old GT (1, 2, 3) to new (reordered) GT (2, 3, 1) for instance.
    gvs_alleles = alt_allele_string.split(",")
    if (len(gvs_alleles) == 0):
        print(f"No ALT ALLELES???")
        sys.exit(1)

    gvs_old_allele_index_to_new = {}
    gvs_old_allele_index_to_new["0"] = "0"
    if (len(gvs_alleles) > 1):
        old_alleles = gvs_alleles[:]
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

def info_diff(gvs_dict, vds_dict, key, locus, gvs_info, vds_info):
    gvs_value = gvs_dict[key]
    vds_value = vds_dict[key]
    if (gvs_value != vds_value):
        print(f"DIFF: Locus {locus}, INFO Field Value for {key} differs between gvs and vds: {gvs_value} vs {vds_value}\nGVS INFO: {gvs_info}\nVDS INFO: {vds_info}")
        if not args.report_all_diffs:
            sys.exit(1)

def simple_diff(gvs_dict, vds_dict, key, locus, sample):
    gvs_value = gvs_dict[key]
    vds_value = vds_dict[key]
    if (gvs_value != vds_value):
        print(f"DIFF: Locus {locus}, Sample {sample} Value for {key} differs between gvs and vds:\n{gvs_value} vs {vds_value}")
        if not args.report_all_diffs:
            sys.exit(1)

def get_sample_genotype_fields(format_fields_to_column, sample_genotype_string):
    tokens = sample_genotype_string.split(":")
    sample_genotype_fields = {}
    for key, index in format_fields_to_column.items():
        value = "."         # Missing field
        if (index < len(tokens)):
            value = tokens[index]
        sample_genotype_fields[key] = value
    return sample_genotype_fields

def get_format_fields_to_column(token):
    tokens = token.split(":")
    fields_to_column = {}
    for i in range(0, len(tokens)):
        fields_to_column[tokens[i]] = i
    return fields_to_column

def get_sample_to_column(fline):
    tokens = fline.split()
    sample_to_column = {}
    for i in range(9, len(tokens)):
        sample_to_column[tokens[i]] = i
    return sample_to_column

def skip_filtered_lines(fline, fp):
    ftokens = fline.split()
    filter = ftokens[6]
    while True:
        if filter == "." or filter == "PASS":
            return fline
        fline = fp.readline().rstrip()
        ftokens = fline.split()
        filter = ftokens[6]


def skip_header(fp):
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

    print(f"NOTE: If a VDS sample genotype record does not contain a FT field (it is missing or explicitly set to '.') then we are considering that a PASS")

    compare_gvs_vcf_with_vds_vcf(args.gvs_vcf, args.vds_vcf, args.skip_gvs_filtered_lines)
