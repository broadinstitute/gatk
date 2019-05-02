import argparse
import logging
import os
import h5py


""" This script copies the hd5 groups specified as 'groups' from all hd5 files within the 'sources'
    directories to the same-named files within the 'destination' directory.
    
    Each source directory in 'sources' must contain the group in 'groups', respectively.
    
    If the destination directory and/or file(s) don't exist, it creates them.
    
    If any of the destination files contains the specified group already, it errors out.
    
    Example command line:
    python .merge_hd5s.py \
        --groups continuous categorical \
        --sources /path/to/src/continuous/tensor/directory /path/to/src/categorical/tensor/directory \
        --dest /path/to/output/directory \
        --logging_level DEBUG
"""

TENSOR_EXT = '.hd5'


def _copy_group(group_name: str, src_file_path: str, dest_file_handle: h5py.File):
    with h5py.File(src_file_path, 'r') as f_src:
        group = f"/{group_name}"

        # Get the name of the parent for the group we want to copy
        group_path = f_src[group].parent.name

        # If the group doesn't exist in the destination, create it (along with parents, if any)
        group_id = dest_file_handle.require_group(group_path)

        # Copy source:/categorical to dest:/categorical
        f_src.copy(group, group_id, name=group)

        msg = f"Copied hd5 group '{group}' from source '{src_file_path}' to '{dest_file_handle.filename}'..."
        logging.debug(msg)


def copy_groups(group_name: str, src_dir: str, dest_dir: str):
    msg_attempting = f"Attempting to copy hd5 files from '{src_dir}' to '{dest_dir}'... "
    logging.debug(msg_attempting)

    if not os.path.exists(src_dir):
        raise ValueError('Source directory does not exist: ', src_dir)

    # If dest_dir doesn't exist, create it
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    # Iterate over files in src_dir
    for src_file in os.listdir(src_dir):
        if src_file.endswith(TENSOR_EXT):
            # Name the destination file with the same as the source's
            src_file_path = os.path.join(src_dir, src_file)
            dest_file_path = os.path.join(dest_dir, src_file)
            with h5py.File(dest_file_path, 'a') as f_dest:
                _copy_group(group_name, src_file_path, f_dest)
        else:
            continue

    msg_succeeded = f"Successfully copied the group '{group_name}' from hd5 files in '{src_dir}' to '{dest_dir}'... "
    logging.debug(msg_succeeded)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--sources', nargs='+', help='List of source directories with hd5 files')
    parser.add_argument('--destination', help='Destination directory to copy hd5 groups to')
    parser.add_argument('--groups', nargs='+')
    parser.add_argument("--logging_level", default='INFO', help="Logging level",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    logging.getLogger().setLevel(args.logging_level)

    for group, source in zip(args.groups, args.sources):
        copy_groups(group, source, args.destination)
