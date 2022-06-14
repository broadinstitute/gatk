import csv
import argparse

expected_subpopulations = [
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
    csvout = csv.writer(csvout, delimiter='\t')
    next(tsvin)  # Skip header row

    for row in tsvin:
      csvout.writerow([row[0], row[4]])


if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract subpopulation per sample data out of a callset TSV')
  parser.add_argument('--input_path',type=str, metavar='path', help='path to the original callset TSV', required=True)
  parser.add_argument('--output_path',type=str, metavar='path', help='path for the output TSV', required=True)

  args = parser.parse_args()

  extract_subpopulation(args.input_path,
                        args.output_path)
