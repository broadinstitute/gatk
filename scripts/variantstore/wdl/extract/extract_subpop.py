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

def extract_subpopulation(input_path):
  seen_subpopulations = []
  with open(input_path, newline='') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    next(tsvin) # skip the header
    for row in tsvin:
      row_subpopulation = row[4]
      if row_subpopulation in expected_subpopulations:
        if row_subpopulation in seen_subpopulations:
          # open the correct output file
          output_file = open(row_subpopulation + "_subpopulation.args", 'a', newline='')
          csvout = csv.writer(output_file, delimiter='\t')
          # add sample name to the subpopulation file
          csvout.writerow([row[0]])
        else:
          seen_subpopulations.append(row_subpopulation)
          # Create a new file for this subpopulation
          output_file = open(row_subpopulation + "_subpopulation.args", 'w', newline='')
          csvout = csv.writer(output_file, delimiter='\t')
          # add sample name to the subpopulation file
          csvout.writerow([row[0]])
          print(seen_subpopulations)
      else:
        print("WARNING: This list had an unexpected subpopulation", row_subpopulation)
      ## TODO do I need to close these?

if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract subpopulation per sample data out of a callset TSV')
  parser.add_argument('--input_path',type=str, metavar='path', help='path to the original callset TSV', required=True)

  args = parser.parse_args()

  extract_subpopulation(args.input_path)

