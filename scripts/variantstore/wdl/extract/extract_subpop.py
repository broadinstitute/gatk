import csv
import argparse

valid_subpopulations = [
 "afr",
 "amr",
 "eas",
 "eur",
 "mid",
 "oth",
 "sas"
]

def extract_subpopulation(input_path, output_path):
  with open(input_path, newline='') as tsvin, open(output_path, 'w', newline='') as csvout:
    tsvin = csv.reader(tsvin, delimiter='\t')
    csvout = csv.writer(csvout, delimiter='\t', lineterminator="\n")
    next(tsvin)  # Skip header row

    observed_subpopulations = set()
    for row in tsvin:
      if row[4] not in valid_subpopulations:
          raise ValueError(f"Unrecognized subpopulation: {row[4]} in {args.input_path}")
      observed_subpopulations.add(row[4])
      csvout.writerow([row[0], row[4]])

    return observed_subpopulations

def write_custom_annotations_files(observed_subpopulations, custom_annotations_template_path):

    with open(custom_annotations_template_path, 'w', newline='') as csvout:
        csvout = csv.writer(csvout, delimiter='\t', lineterminator="\n")
        csvout.writerow(["#title=gvsAnnotations"])
        csvout.writerow(["#assembly=GRCh38"])
        csvout.writerow(["#matchVariantsBy=allele"])

        chrom_line = ["#CHROM", "POS", "REF", "ALT", "AC", "AN", "AF", "AC_Hom", "AC_Het", "AC_Hemi"]
        categories_line = ["#categories", ".", ".", ".", "AlleleCount", "AlleleNumber", "AlleleFrequency", "AlleleCount", "AlleleCount", "AlleleCount"]
        description_line = ["#descriptions", ".", ".", ".", ".", ".", ".", ".", ".", "."]
        type_line = ["#type", ".", ".", ".", "number", "number", "number", "number", "number", "number"]

        for subpopulation in sorted(observed_subpopulations):
            for annotation in (["AC", "AN", "AF", "AC_Hom", "AC_Het", "AC_Hemi"]):
                chrom_line.append(f"{annotation}_{subpopulation}")
            categories_line += ["AlleleCount", "AlleleNumber", "AlleleFrequency", "AlleleCount", "AlleleCount", "AlleleCount"]
            for i in range(6):
                description_line.append(".")
                type_line.append("number")

        csvout.writerow(chrom_line)
        csvout.writerow(categories_line)
        csvout.writerow(description_line)
        csvout.writerow(type_line)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract subpopulation per sample data out of a callset TSV')
  parser.add_argument('--input_path', type=str, metavar='path', help='path to the original callset TSV', required=True)
  parser.add_argument('--output_path', type=str, metavar='path', help='path for the output TSV', required=True)
  parser.add_argument('--custom_annotations_template_path', type=str, metavar='path', help='path for the custom annotations template file', required=True)

  args = parser.parse_args()

  observed_subpopulations = extract_subpopulation(args.input_path,
                                                  args.output_path)
  write_custom_annotations_files(observed_subpopulations,
                                 args.custom_annotations_template_path)
