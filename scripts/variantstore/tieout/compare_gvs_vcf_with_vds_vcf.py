import argparse
import sys

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

                if (gvs_tokens[2] != vds_tokens[2]):
                    print(f"DIFF: ID differs between VCF line:\n{gvs_tokens[2]} vs {vds_tokens[2]}")
                    if not args.report_all_diffs:
                        sys.exit(1)

                if (gvs_tokens[3] != vds_tokens[3]):
                    print(f"DIFF: REF differs between VCF line:\n{gvs_tokens[3]} vs {vds_tokens[3]}")
                    if not args.report_all_diffs:
                        sys.exit(1)

                # Alt Alleles!
                gvs_alleles = gvs_tokens[4].split(",")
                if (len(gvs_alleles) > 1):
                    gvs_alleles.sort()
                    gvs_tokens[4] = ",".join(gvs_alleles)

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

                # TODO - Skipping INFO for now
                # if (gvs_tokens[7] != vds_tokens[7]):
                #     print(f"DIFF: INFO differs between VCF line:\n{gvs_tokens[7]} vs {vds_tokens[7]}")
                #     sys.exit(1)

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

                    # NOTE: More complex rules for GT
                    # simple_diff(gvs_sample_genotypes, vds_sample_genotypes, "GT", locus, sample)
                    gvs_gt = gvs_sample_genotypes["GT"]
                    vds_gt = vds_sample_genotypes["GT"]
                    vds_lgt = vds_sample_genotypes["LGT"]
                    if (vds_gt == "./." and vds_ft == "FAIL"):
                        # NOTE/HMM. In the VDS we no call a genotype if the FT is a FAIL, the original genotype is in LGT.
                        vds_gt = vds_lgt
                    if (gvs_gt != vds_gt):
                        print(f"\nDIFF: Locus {locus}, Sample {sample} Value for GT differs between gvs and vds: {gvs_gt} vs {vds_gt}")
                        print(f"Locus {locus}.  gvs FORMAT: {gvs_tokens[8]} vds FORMAT: {vds_tokens[8]}")
                        print(f"Sample {sample}.    gvs sample: {gvs_tokens[gvs_sample_to_column[sample]]}  vds_sample: {vds_tokens[vds_sample_to_column[sample]]}")
                        if not args.report_all_diffs:
                            sys.exit(1)

                    simple_diff(gvs_sample_genotypes, vds_sample_genotypes, "GQ", locus, sample)
                    simple_diff(gvs_sample_genotypes, vds_sample_genotypes, "RGQ", locus, sample)

                    ## AD
                    # gvs_ad = gvs_sample_tokens[gvs_format_fields_to_column["AD"]]
                    # vds_ad = vds_sample_tokens[vds_format_fields_to_column["LAD"]]
                    gvs_ad = gvs_sample_genotypes["AD"]
                    vds_ad = vds_sample_genotypes["LAD"]
                    if (gvs_ad != vds_ad):
                        print(f"\nDIFF: Locus {locus}, Sample {sample} Value for differs between gvs key AD and vds key LAD: {gvs_ad} vs {vds_ad}")
                        print(f"Locus {locus}.  gvs FORMAT: {gvs_tokens[8]} vds FORMAT: {vds_tokens[8]}")
                        print(f"Sample {sample}.    gvs sample: {gvs_tokens[gvs_sample_to_column[sample]]}  vds_sample: {vds_tokens[vds_sample_to_column[sample]]}")
                        if not args.report_all_diffs:
                            sys.exit(1)

                    # if (gvs_sample != vds_sample):
                    #     print(f"DIFF: SAMPLE differs between VCF line:\n{gvs_sample} vs {vds_sample}")
                        # sys.exit(1)

                gvs_line = gvs.readline().rstrip()
                vds_line = vds.readline().rstrip()
                # if count > 10000:
                #     sys.exit(1)

            exit(0)

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

    print(f"NOTE/TODO: If a VDS sample genotype record does not contain a FT field (it is missing or explicitly set to '.') then we are considering that a PASS")
    print(f"NOTE/TODO: Not doing ANYTHING with the gds VCF's sample genotype record 'LA' field\n")

    compare_gvs_vcf_with_vds_vcf(args.gvs_vcf, args.vds_vcf, args.skip_gvs_filtered_lines)
