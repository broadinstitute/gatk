# ML4H is released under the following BSD 3-Clause License:
#
# Copyright (c) 2020, Broad Institute, Inc. and The General Hospital Corporation.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name Broad Institute, Inc. or The General Hospital Corporation
#   nor the names of its contributors may be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
import time
import json
from multiprocessing import Pool, cpu_count
from functools import partial
import os
import numpy as np
import pandas as pd
from typing import List, Tuple, Optional, Dict
from collections import defaultdict
import xml.etree.ElementTree as ET
import string
import uuid

# Set of printable characters.
printable = set(string.printable)

# Constant: ignore these keys in the XML during ingestion.
# For the UK Biobank ECG XMLs:
ukb_ignore_elements = [
    "Stripdata",
    "Fulldisclosure",
    "Mediansamples",
    "Waveformdata",
    "Fulldisclosuredata",
]
# For Partners (MUSE) ECG XMLs:
muse_ignore_elements = ['MeasurementMatrix','WaveFormData']

# For the UKB ECG XMLs
def _sample_id_from_path(path):
    # Assumes thef format UKBID_FIELDID_INSTANCE.csv
    base = os.path.split(path)[-1].split(".")[0].split("_")
    return base[0], base[2]  # ID, INSTANCE


# For MUSE ECG XMLs
# def _sample_id_from_path(path):
#     # base = ?
#     return '%032x' % random.randrange(16**32),0 # random, nonsense


def generate_output_name(
    subject_id: str, instance: str, source: str = "ukb", data_type: str = "ecg"
) -> str:
    """Helper function to generate output names. Current implementation will
    return `BROAD_ML4H_mdrk_{source}_{data_type}_metadata_{id}_i{instance}_{hash}`.

    Examples:
        * BROAD_ML4H_mdrk_ukb_ecg_metadata_6025769_i2_cc96e96cbe03272703223adf
        * BROAD_ML4H_mdrk_muse_ecg_metadata_2331765_i0_545d840638e5fe00592ee64d.pq

    Args:
        subject_id (str): MUSE or UKB subject identifier. This could also be any random value for privacy.
        instance (str): UKB instance number. For MUSE this can be anything.
        source (str, optional): XML source: currently `ukb` or `muse` but no restrictions apply. Defaults to 'ukb'.
        data_type (str, optional): Modality associated with the XML. Defaults to 'ecg'.

    Returns:
        str: A prepared output string name.
    """
    return f"BROAD_ML4H_mdrk_{str(source)}_{data_type}_metadata_{str(subject_id)}_i{str(instance)}_{uuid.uuid4().hex}"


def _recurseTree(root, path: str, store_dict: dict, ignore_elements: list):
    """Recursive function for iterating over the XML tree and extracting out all
    keys and values for those keys that are not in the `ignore_elements` black list.
    Note: All values will be stored as strings and all strings are scrubbed for
    newlines, tabs, and trailing and leading white space. Multiple elements for the
    same putative path will be stored as an array of the values for the same key.
    Nones are stored as empty strings ('').

    Args:
        root: Current element in the XML tree.
        path (str): Current string name for the traversed path.
        store_dict (dict): Output dictionary.
        ignore_elements (list): List of keys to ignore
    """
    if root.tag in ignore_elements or root.tag.title() in ignore_elements:
        return

    data = root.attrib.get("name", root.text)
    if data is not None:  # Not None
        data = (
            "".join(filter(lambda x: x in printable, data))
            .replace("\n", "")
            .replace("\r", "")
            .replace("\t", "")
            .strip()
        )

        key = path + "_" + root.tag.title()
        if key in store_dict:  # Key exists
            if isinstance(store_dict[key], list) == False:
                store_dict[key] = [store_dict[key]]
            store_dict[key] = store_dict[key] + [data]
        else:  # Key does not exist
            store_dict[key] = data
    else:  # Is none
        key = path + "_" + root.tag.title()
        if key in store_dict:  # Key exist
            if isinstance(store_dict[key], list) == False:
                store_dict[key] = [store_dict[key]]
            store_dict[key] = store_dict[key] + [""]
        else:  # Key does not exist
            store_dict[key] = ""

    # Recurse over children.
    for elem in root.getchildren():
        _recurseTree(
            elem, path + "_" + root.tag.title(), store_dict, ignore_elements
        )


def ingest_metadata_from_xml(
    file: str,
    sample_id: str,
    instance: str,
    destination: str,
    ignore_elements: list,
    data_source: str = "ukb",
    data_type: str = "ecg",
):
    """Extract all elements (keys) and their values from an XML through recursion. Keys
    in the global `ukb_ignore_elements` list will be ignored and not extracted.

    Args:
        file (str): [description]
        sample_id (str): [description]
        instance (str): [description]
        destination (str): [description]
        ignore_elements (list):
        data_source (str):
        data_type (str, optional): [description]. Defaults to 'ecg'.
    """
    tree = ET.ElementTree(file=file)  # XML tree
    root = tree.getroot()  # XML root
    path = ""  # Starting path
    store_dicter = {}  # Output dictionary
    _recurseTree(root, path, store_dicter, ignore_elements)

    # Lists of values must be stored as lists of lists to be compatible
    # with Pandas.
    for k in store_dicter.keys():
        if isinstance(store_dicter[k], list):
            store_dicter[k] = [store_dicter[k]]

    # Convert dictionary into a Pandas DataFrame and cast all columns as strings.
    # This will also cast lists as string representation of strings.
    df = pd.DataFrame(store_dicter).astype(str)
    # Keep source file path for bookkeeping.
    df["source_file"] = file
    df.to_parquet(
        os.path.join(
            destination,
            generate_output_name(
                sample_id, instance, source=data_source, data_type=data_type
            )
            + ".pq",
        ),
        compression="zstd",
    )


def _process_file(
    path: str,
    destination: str,
    ignore_elements: list,
    data_source: str = "ukb",
    data_type: str = "ecg",
) -> Tuple[str, Optional[str]]:
    sample_id, instance = _sample_id_from_path(path)
    try:
        ingest_metadata_from_xml(
            path,
            sample_id,
            instance,
            destination,
            ignore_elements,
            data_source,
            data_type,
        )
        return path, None
    except Exception as e:
        return path, str(e)


def _process_files(
    files: List[str],
    destination: str,
    ignore_elements: list,
    data_source: str = "ukb",
    data_type: str = "ecg",
) -> Dict[str, str]:
    errors = {}
    name, _ = _sample_id_from_path(files[0])
    process_file = partial(
        _process_file,
        destination=destination,
        data_source=data_source,
        data_type=data_type,
        ignore_elements=ignore_elements,
    )

    print(f"Starting process {name} with {len(files)} files")
    for i, (path, error) in enumerate(map(process_file, files)):
        if error is not None:
            errors[path] = error

        if len(files) % max(i // 10, 1) == 0:
            print(f"{name}: {(i + 1) / len(files):.2%} done")

    return errors


def _partition_files(files: List[str], num_partitions: int) -> List[List[str]]:
    """Split files into num_partitions partitions of close to equal size"""
    id_to_file = defaultdict(list)
    for f in files:
        id_to_file[_sample_id_from_path(f)[0]].append(f)
    sample_ids = np.array(list(id_to_file))
    np.random.shuffle(sample_ids)
    split_ids = np.array_split(sample_ids, num_partitions)
    splits = [
        sum((id_to_file[sample_id] for sample_id in split), []) for split in split_ids
    ]
    return [split for split in splits if split]  # lose empty splits


def multiprocess_ingest(
    files: List[str],
    destination: str,
    data_source: str = "ukb",
    data_type: str = "ecg",
    ignore_elements: list = ukb_ignore_elements,
):
    """Embarassingly parallel ingestion wrapper.

    Example usage:
    ```
    >>> import glob
    >>> files = glob.glob('/some/path/xmls/*.xml')
    >>> multiprocess_ingest(files, '/output_path/level1')
    ```

    Args:
        files (List[str]): Input list of files.
        destination (str): Output destination on disk.
        data_source (str): Data source (e.g. 'ukb'). Defaults to 'ukb'
        data_type (str): Data source type (e.g. 'ecg'). Defaults to 'ecg'

    Returns:
        [dict]: Returns a dictionary of encountered errors.
    """
    print(f"Beginning ingestion of {len(files)} files using {cpu_count()} threads.")
    os.makedirs(destination, exist_ok=True)
    # Partition files to prevent races.
    split_files = _partition_files(files, cpu_count())

    errors = {}
    start = time.time()
    with Pool(cpu_count()) as pool:
        results = [
            pool.apply_async(
                _process_files,
                (split, destination, ignore_elements, data_source, data_type),
            )
            for split in split_files
        ]
        for result in results:
            errors.update(result.get())

    delta = time.time() - start
    print(f"Ingestion took {delta:.1f} seconds at {delta / len(files):.1f} s/file")

    with open(os.path.join(destination, "errors.json"), "w") as f:
        json.dump(errors, f)

    return errors
