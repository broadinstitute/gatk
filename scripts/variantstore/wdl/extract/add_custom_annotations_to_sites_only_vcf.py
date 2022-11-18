import argparse
import sys

def add_custom_annotations_to_sites_only_vcf(sites_only_vcf, annot_tsv, output_vcf):
    with open(annot_tsv) as tsv:
        names, names_to_column = skip_annot_header(tsv, annot_tsv)
        with open(sites_only_vcf) as vcf:
            vcf_header, chrom_line = skip_vcf_header(vcf)
            if not chrom_line.endswith("INFO"):
                print(f"ERROR: Sites only VCF {sites_only_vcf} has unexpected header")
                sys.exit(1)

            with open(output_vcf, 'w') as out:
                out.write(vcf_header)

                count = 0
                while True:
                    count += 1

                    vcf_line, annot_line = co_iterate(vcf, tsv)
                    tokens = vcf_line.split("\t")

                    names_to_value = get_annotation_map(annot_line, names)
                    info_fields = []
                    for key, value in names_to_value.items():
                        if (key != "#CHROM" and key != "POS" and key != "REF" and key != "ALT"):
                            token = key + "=" + value
                            info_fields.append(token)
                    info_field = ';'.join(info_fields)

                    if (tokens[0] != names_to_value["#CHROM"] or tokens[1] != names_to_value["POS"] or tokens[3] != names_to_value["REF"] or tokens[4] != names_to_value["ALT"]):
                        print(f"Mismatch between keys in VCF line:\n{vcf_line} and those in annotation line:\n{annot_line}")
                        sys.exit(1)

                    if (tokens[7] != "."):
                        # TODO - we really should be able to handle this easily, but don't expect to see this situation.
                        print(f"ERROR: Sites only VCF {sites_only_vcf} has values in INFO field")
                        sys.exit(1)

                    tokens[7] = info_field
                    vcf_line = '\t'.join(tokens) + '\n'
                    out.write(vcf_line)

def co_iterate(vcf_fp, tsv_fp):
    annot_line = tsv_fp.readline().rstrip()
    vcf_line = vcf_fp.readline().rstrip()
    if (annot_line == "" and vcf_line == ""):
        print("Done")
        sys.exit(0)
    tokens = vcf_line.split("\t")
    if tokens[4] == "*":
        print(f"WARNING: Skipping Spanning deletion on VCF line:\n{vcf_line}")
        vcf_line = vcf_fp.readline().rstrip()
    return(vcf_line, annot_line)



def skip_vcf_header(fp):
    """
    A method to skip over the header in the VCF
    :param fp: The file's file pointer
    :return: The nextmost line that is not a header line
    """
    vcf_header = ""
    while (line := fp.readline().rstrip()):
        vcf_header += line + '\n'
        if not line.startswith("##"):
            return vcf_header, line

def get_annotation_map(line, names_array):
    tokens = line.split("\t")
    if (len(tokens) != len(names_array)):
        print(f"ERROR: Different number of values in line: {line} than found in header")
    names_to_values = {}
    for i in range(0, len(tokens)):
        names_to_values[names_array[i]] = tokens[i]
    return names_to_values

def skip_annot_header(fp, annot_tsv):
    """
    A method to skip over the header in the custom annotations file
    AND returns a map of the column names to their index
    :param fp: The file's file pointer
    :return: a dictionary of the field names to their index in the line
    """

    fp.readline()
    fp.readline()
    fp.readline()
    line = fp.readline().rstrip()
    if not line.startswith("#CHROM"):
        print(f"ERROR: Fourth line of {annot_tsv} does NOT start with '#CHROM!'")
        sys.exit(1)
    names_to_column = get_name_to_column_map(line)
    # Skip over the other header lines.
    fp.readline()
    fp.readline()
    fp.readline()

    return names_to_column


def get_name_to_column_map(line):
    """
    A method to parse the #CHROM line of the custom annotations TSV file and return a map of
    name to column index.
    :param line: The #CHROM line from the custom annotations TSV file header
    :return: a dictionary of the field names to their index in the line
    """
    tokens = line.split("\t")
    names_to_column = {}
    for i in range(0, len(tokens)):
        names_to_column[tokens[i]] = i
    return tokens, names_to_column


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='A tool to parse a custom annotations file and add those annotations as INFO fields to a sites-only VCF')

    parser.add_argument('--input_vcf', type=str, help='Input (sites-only) VCF', required=True)
    parser.add_argument('--custom_annotations_tsv', type=str, help='Custom annotations file', required=True)
    parser.add_argument('--output_vcf', type=str, help='Output (sites-only) VCF', required=True)
    parser.add_argument("--verbose", action="store_true", help="increase output verbosity")
    args = parser.parse_args()

    add_custom_annotations_to_sites_only_vcf(args.input_vcf, args.custom_annotations_tsv, args.output_vcf)