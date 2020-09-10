import argparse
import os
import csv
from getpass import getuser
from collections import defaultdict

JOIN_CHAR = '_'
SCRIPT_PATH = 'tensormap/mgb/tensor_maps_partners_ecg_labels.py'
TENSOR_FROM_FILE_MAKER = "make_partners_ecg_label"
PREFIX = "partners_ecg"
TENSOR_PATH_PREFIX = "partners_ecg_rest"


def _clean_label_string(string):
    """Replace spaces and slashes with JOIN_CHAR,
       and remove parentheses and commas"""
    string = string.replace(' ', JOIN_CHAR)
    string = string.replace('/', JOIN_CHAR)
    string = string.replace('(', '')
    string = string.replace(')', '')
    string = string.replace(',', '')
    return string


def _write_tmap_to_py(write_imports, py_file, label_maps, channel_maps, keys_in_hd5):
    """Given label_maps (which associates labels with source phrases)
       and channel_maps (which associates labels with unique sublabels),
       define the tensormaps to associate source phrases with precise labels,
       and write these maps in a python script
    """

    # Add import statements to .py
    if write_imports:
        py_file.write(f"from ml4h.TensorMap import TensorMap, Interpretation\n")
        py_file.write(f"from ml4h.tensormap.mgb.ecg import {TENSOR_FROM_FILE_MAKER}\n\n")

    for label in label_maps:
        cm = '{'

        for i, channel in enumerate(channel_maps[label]):
            cm += f"'{channel}': {i}, "

        # At this point, i = len(channel_maps[label])-1
        # If 'unspecified' is not a label, we need to add and index it
        if 'unspecified' not in channel_maps[label]:
            cm += f"'unspecified': {i+1}"

        cm += '}'

        for key in keys_in_hd5:
            py_file.write(f"{PREFIX}_{key}_{label} = TensorMap('{PREFIX}_{key}_{label}', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='{TENSOR_PATH_PREFIX}', channel_map={cm}, tensor_from_file={TENSOR_FROM_FILE_MAKER}(keys='{key}', dict_of_list = {label_maps[label]})) \n\n")
            py_file.write(f"{PREFIX}_{key}_{label}_newest = TensorMap('{PREFIX}_{key}_{label}_newest', interpretation=Interpretation.CATEGORICAL, path_prefix='{TENSOR_PATH_PREFIX}', channel_map={cm}, tensor_from_file={TENSOR_FROM_FILE_MAKER}(keys='{key}', dict_of_list = {label_maps[label]})) \n\n")


def _write_partners_ecg_tmap_script(py_file, partners_ecg_label_dir, keys_in_hd5):
    # Set flag for writing import statements to true
    write_imports = True

    # Iterate through all files in the mgb CSV labels folder
    for file in os.listdir(partners_ecg_label_dir):

        # Ignore files that do not end in .csv
        if not file.endswith('.csv'):
            continue

        # Isolate the task name
        task = file.replace('.csv', '').replace('c_', '')

        # Create list of lists;
        # outer list is rows in CSV,
        # inner list is columns; idx 0 is source phrase
        #                        idx 1 is label (level 1 of hierarchy)
        #                        idx 2 is label (level 2 of hierarchy) etc.
        fpath_this_task = os.path.join(partners_ecg_label_dir, file)
        lol = list(csv.reader(open(fpath_this_task), delimiter=','))

        # Associate labels with source phrases in dict of dicts:
        #   keys   - task name and all oot-level labels in hierarchy
        #   values - dicts:
        #       keys   - labels (next level down in hierarchy)
        #       values - list of source phrases that map to a given label
        # Note: because the first key is the task name, keys of dicts in
        # label_map[task] are the remaining keys in label_map itself
        label_maps = defaultdict(dict)

        # Associate labels with unique set of sublabels in dict of sets
        # keys   - every label in hierarchy with children
        # values - set of all child labels within a given label
        channel_maps = defaultdict(set)

        # Iterate through every source phrase in list of lists (label map)
        for row in lol:

            # Initialize prefix as empty list
            prefix = []

            # For the row, iterate through label strings 2:end
            for label_str in row[1:]:

                # If the label string is blank, skip
                if label_str == '':
                    continue

                # Clean label string
                label_str = _clean_label_string(label_str)

                # If ???
                if len(prefix) == 0:
                    channel_maps[task].add(label_str)

                    # If we already have a list of source phrases for this
                    # task and label string, append the source phrase to it
                    if label_str in label_maps[task]:
                        label_maps[task][label_str].append(row[0])

                    # If no list exists yet for this task and label string,
                    # create a new list and initialize with the source phrase
                    else:
                        label_maps[task][label_str] = [row[0]]
                else:
                    # Join elements in prefix list by JOIN_CHAR
                    prefix_merged = JOIN_CHAR.join(prefix)
                    channel_maps[prefix_merged].add(label_str)

                    if label_str in label_maps[prefix_merged]:
                        label_maps[prefix_merged][label_str].append(row[0])
                    else:
                        label_maps[prefix_merged][label_str] = [row[0]]

                    label_maps[prefix_merged][label_str].append(row[0])

                prefix.append(label_str)

        _write_tmap_to_py(write_imports, py_file, label_maps, channel_maps, keys_in_hd5)
        write_imports = False

        print(f"Created Tensor Maps from {file} and saved in {SCRIPT_PATH}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--label_maps_dir', default=f'/home/{getuser()}/partners_ecg_labeling', help='Path to directory with c_task.csv label maps.')
    parser.add_argument('--hd5_keys', default=['read_md_clean', 'read_pc_clean'], nargs='+')
    args = parser.parse_args()

    # Open the python script, parse the .csv label maps,
    # and create code that creates TensorMaps
    with open(SCRIPT_PATH, 'w') as py_file:
        _write_partners_ecg_tmap_script(py_file, args.label_maps_dir, args.hd5_keys)
