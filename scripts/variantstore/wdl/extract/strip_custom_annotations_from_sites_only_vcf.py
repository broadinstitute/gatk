import argparse
import sys

def strip_custom_annotations_from_sites_only_vcf(input_sites_only_vcf, input_custom_annotations_tsv, output_vcf, output_custom_annotations_tsv):
    with open(input_custom_annotations_tsv) as tsv:
        header, column_names, names_to_column_map = parse_annot_header(tsv, input_custom_annotations_tsv)
        with open(output_custom_annotations_tsv, 'w') as out_tsv:
            out_tsv.write(header)
            with open(input_sites_only_vcf) as input_vcf:
                vcf_header, chrom_line = skip_vcf_header(input_vcf)
                with open(output_vcf, 'w') as out_vcf:
                    out_vcf.write(vcf_header)

                    if not chrom_line.endswith("INFO"):
                        print(f"ERROR: Input Sites only VCF {input_sites_only_vcf} has unexpected header")
                        sys.exit(1)

                    while True:
                        vcf_line = input_vcf.readline().rstrip()
                        if (vcf_line == ""):
                            print("Done")
                            sys.exit(0)

                        vcf_fields = vcf_line.split("\t")
                        info_field_map = get_info_fields(vcf_fields[7])

                        # Populate the annotations line
                        annotations = []
                        annotations.append(vcf_fields[0])   # Chrom
                        annotations.append(vcf_fields[1])   # Position
                        annotations.append(vcf_fields[3])   # Ref
                        annotations.append(vcf_fields[4])   # Alt
                        for index in range(0, len(column_names)):
                            key = column_names[index]
                            if (key in info_field_map):
                                value = info_field_map[key]
                                del info_field_map[key]
                            else:
                                value = "."
                            annotations.append(value)

                        # Error checking - did we find any unexpected fields in the VCF INFO field?
                        if (len(info_field_map) > 0):
                            print(f"ERROR: Found fields in the VCF INFO field (not found in annotations file)")
                            for key, value in info_field_map.items():
                                print(f"{key}={value}")
                            sys.exit(1)

                        out_vcf.write("\t".join(vcf_fields[0:-1]) + "\t.\n")
                        out_tsv.write("\t".join(annotations) + "\n")

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


def skip_vcf_header(fp):
    """
    A method to skip over the header in the VCF
    :param fp: The file's file pointer
    :return: The VCF Header (as a text block) and th next line that is not a header line
    """
    vcf_header = ""
    for line in fp:
        line = line.rstrip()
        vcf_header += line + '\n'
        if not line.startswith("##"):
            return vcf_header, line


def parse_annot_header(fp, annot_tsv):
    """
    A method to skip over the header in the custom annotations file
    AND returns a map of the column names to their index
    :param fp: The file's file pointer
    :return: a dictionary of the field names to their index in the line
    """

    annot_header = ""
    annot_header += fp.readline().rstrip() + '\n'
    annot_header += fp.readline().rstrip() + '\n'
    annot_header += fp.readline().rstrip() + '\n'
    line = fp.readline().rstrip()
    annot_header += line + '\n'
    if not line.startswith("#CHROM"):
        print(f"ERROR: Fourth line of {annot_tsv} does NOT start with '#CHROM!'")
        sys.exit(1)
    names, names_to_column = get_name_to_column_map(line)
    if (names[0] != "#CHROM" or names[1] != "POS" or names[2] != "REF" or names[3] != "ALT"):
        print(f"Unexpected column names in {annot_tsv}")
        sys.exit(1)
    names = names[4:]

    # Skip over the other header lines.
    annot_header += fp.readline().rstrip() + '\n'
    annot_header += fp.readline().rstrip() + '\n'
    annot_header += fp.readline().rstrip() + '\n'

    return annot_header, names, names_to_column


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
    parser = argparse.ArgumentParser(allow_abbrev=False, description='A tool to take the custom annotations added to a sites-only VCF and recreate a custom annotations file from them')

    parser.add_argument('--input_vcf', type=str, help='Input (sites-only) VCF', required=True)
    parser.add_argument('--input_custom_annotations_tsv', type=str, help='Input custom annotations header file', required=True)
    parser.add_argument('--output_vcf', type=str, help='Output (sites-only) VCF - will be stripped of custom annotations', required=True)
    parser.add_argument('--output_custom_annotations_tsv', type=str, help='Output custom annotations file', required=True)
    parser.add_argument("--verbose", action="store_true", help="increase output verbosity")
    args = parser.parse_args()

    strip_custom_annotations_from_sites_only_vcf(args.input_vcf, args.input_custom_annotations_tsv, args.output_vcf, args.output_custom_annotations_tsv)