import logging
import os
import h5py

# TODO: Use cmdline args for file paths and log level
CATEGORICAL_TENSORS_SRC_DIR = '/Users/kyuksel/ml4cvd/tensors/dataflow_tensors/tensors_ukbb_dev_categorical'
CONTINUOUS_TENSORS_SRC_DIR = '/Users/kyuksel/ml4cvd/tensors/dataflow_tensors/tensors_ukbb_dev_continuous'
DEST_DIR = '/Users/kyuksel/ml4cvd/tensors/dataflow_tensors/tensors_ukbb_dev_categorical_continuous'

CATEGORICAL_GROUP_NAME = '/categorical'
CONTINUOUS_GROUP_NAME = '/continuous'

LOG_LEVEL = logging.DEBUG
TENSOR_EXT = '.hd5'


def _copy_group(group_name: str, src_file_path: str, dest_file_handle: h5py.File):
    with h5py.File(src_file_path, 'r') as f_src:
        # Get the name of the parent for the group we want to copy
        group_path = f_src[group_name].parent.name

        # If the group doesn't exist in the destination, create it (along with parents, if any)
        group_id = dest_file_handle.require_group(group_path)

        # Copy source:/categorical to dest:/categorical
        f_src.copy(group_name, group_id, name=group_name)

        msg = "Copied hd5 group {} from source {} to {}... ".format(group_name, src_file_path, dest_file_handle.filename)
        logging.debug(msg)


def copy_groups(group_name: str, src_dir: str, dest_dir: str):
    msg_attempting = "Attempting to copy hd5 files from {} to {}... ".format(src_dir, dest_dir)
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

    msg_succeeded = "Successfully copied the hd5 files from {} to {}... ".format(src_dir, dest_dir)
    logging.debug(msg_succeeded)


if __name__ == "__main__":
    logging.getLogger().setLevel(LOG_LEVEL)

    copy_groups(CATEGORICAL_GROUP_NAME, CATEGORICAL_TENSORS_SRC_DIR, DEST_DIR)
    copy_groups(CONTINUOUS_GROUP_NAME, CONTINUOUS_TENSORS_SRC_DIR, DEST_DIR)
