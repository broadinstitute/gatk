#!/bin/python

# This script aggregates score tables and metrics files from multiple PathSeq runs
# into a small set of tabular text files more amenable to comparative analysis.

import argparse
import csv

#################
# CONSTANTS
#################

tax_id_column = 'tax_id'
expected_info_columns = [tax_id_column, 'taxonomy', 'type', 'name', 'kingdom', 'reference_length']
expected_score_columns = ['score', 'score_normalized', 'reads', 'unambiguous']

metrics_header_line_symbol = '#'
expected_score_metrics_columns = ['MAPPED_READS', 'UNMAPPED_READS']
expected_filter_metrics_columns = ['PRIMARY_READS', 'READS_AFTER_PREALIGNED_HOST_FILTER', 'READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER', 'READS_AFTER_HOST_FILTER', 'READS_AFTER_DEDUPLICATION', 'FINAL_PAIRED_READS', 'FINAL_UNPAIRED_READS', 'FINAL_TOTAL_READS', 'LOW_QUALITY_OR_LOW_COMPLEXITY_READS_FILTERED', 'HOST_READS_FILTERED', 'DUPLICATE_READS_FILTERED']

##################################
# FUNCTION DEFINITIONS
##################################

# Gets index of column in list header with appropriate exception if it doesn't exist
def get_column_index(header, column_header, file_name):
  try:
    index = header.index(column_header)
    return index
  except:
    raise ValueError("Expected to find column \"" + column_header + "\" in file " + file_name)

# Reads list of score files and returns the score data and taxon information
def read_score_files(file_paths):
  data = {};
  tax_info = {};
  for file_path in file_paths:
    with open(file_path,'rb') as f:
      reader = csv.reader(f, delimiter='\t', strict=True)
      header = reader.next()
      tax_id_index = get_column_index(header, tax_id_column, file_path)
      score_column_indices = []
      for column in expected_score_columns:
        score_column_indices.append(get_column_index(header, column, file_path))
      info_column_indices = []
      for column in expected_info_columns:
        info_column_indices.append(get_column_index(header, column, file_path))
      num_cols = len(header)
      file_data = {}
      for row in reader:
        if len(row) != num_cols:
          raise ValueError("Malformed scores file. Header has " + str(num_cols) + " columns but found a row with " + str(len(row)) + ": " + str(row))
        tax_id = row[tax_id_index]
        if tax_id not in tax_info:
          tax_info[tax_id] = []
          for i in range(len(info_column_indices)):
            tax_info[tax_id].append(row[info_column_indices[i]])
        file_data[tax_id] = []
        for column_index in score_column_indices:
          value = row[column_index]
          file_data[tax_id].append(value)
      data[file_path] = file_data
  return data, tax_info

# Writes aggregated score data for column at data_index in expected_score_columns
def write_score_file(file_paths, base_path, data_index, data, info):
  with open(base_path + '.' + expected_score_columns[data_index] + '.txt', 'w') as f:
    tax_ids = data.keys()
    f.write('\t'.join(expected_info_columns))
    f.write('\t' + '\t'.join(file_paths) + '\n')
    for tax_id in info:
      tax_id_info = info[tax_id]
      f.write('\t'.join(tax_id_info) + '\t')
      for sample in file_paths:
        if tax_id in data[sample]:
          f.write(data[sample][tax_id][data_index] + '\t')
        else:
          f.write('0\t')
      f.write('\n')

# Reads Picard-style metrics file
def read_metrics_file(file_paths, expected_columns):
  metrics = []
  for file_path in file_paths:
    with open(file_path,'rb') as f:
      line = f.readline().strip()
      while (not line) or line[0] == metrics_header_line_symbol:
        line = f.readline().strip()
      header = line.split('\t')
      column_indices = []
      for column in expected_columns:
        column_indices.append(get_column_index(header, column, file_path))
      values = f.readline().strip().split('\t')
      metrics.append([values[i] for i in column_indices])
  return metrics

# Writes aggregated metrics file
def write_metrics_file(file_paths, filter_metrics, score_metrics, base_path):
  with open(base_path + '.metrics.txt', 'w') as f:
    f.write('\t' + '\t'.join(file_paths) + '\n')
    write_metrics_type(f, filter_metrics, expected_filter_metrics_columns)
    write_metrics_type(f, score_metrics, expected_score_metrics_columns)

# Helper for generating table rows
def write_metrics_type(file, metrics, columns):
  for i in range(len(columns)):
    m = [x[i] for x in metrics]
    file.write(columns[i] + '\t' + '\t'.join(m) + '\n')

# Reads list of files from file into an array
def read_files_list(list_path):
  with open(list_path, 'r') as f:
    lines = f.read().splitlines()
    if not lines:
      raise ValueError("List " + list_path + " was empty")
    return lines

####################################
# COMMAND-LINE ARGUMENT DEFINITIONS
####################################

parser = argparse.ArgumentParser(description='Combines PathSeq score and metrics files across multiple samples')
parser.add_argument('output_base_path', help="Base path for output files, e.g. /output_dir/results will generate files /output_dir/results.score.txt, etc.")
parser.add_argument('score_files_list', help="File listing score file paths (1 per line)")
parser.add_argument('filter_metrics_files_list', help="File listing filter metric file paths (1 per line)")
parser.add_argument('score_metrics_files_list', help="File listing score metric file paths (1 per line). Must be the same length as filter_metrics_files_list.")
args = parser.parse_args()

##################################
# MAIN ROUTINE
##################################

score_files = read_files_list(args.score_files_list)
filter_metrics_files = read_files_list(args.filter_metrics_files_list)
score_metrics_files = read_files_list(args.score_metrics_files_list)

if len(filter_metrics_files) != len(score_metrics_files):
  raise ValueError("Lengths of filter_metrics_files_list and score_metrics_files_list must be equal")

score_data, score_info = read_score_files(score_files)
filter_metrics = read_metrics_file(filter_metrics_files, expected_filter_metrics_columns)
score_metrics = read_metrics_file(score_metrics_files, expected_score_metrics_columns)

for i in range(len(expected_score_columns)):
  write_score_file(score_files, args.output_base_path, i, score_data, score_info)

write_metrics_file(filter_metrics_files, filter_metrics, score_metrics, args.output_base_path)
