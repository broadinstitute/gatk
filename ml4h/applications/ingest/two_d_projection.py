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
from typing import Dict, List, Tuple
import os
import h5py
import numpy as np
import pandas as pd
from scipy.ndimage import zoom
from fastparquet import ParquetFile
import matplotlib.pyplot as plt
from ingest_mri import compress_and_store, read_compressed


def project_coronal(x: np.ndarray) -> np.ndarray:
    """Computes the mean 2D projection in the putative coronal dimension
    given axial input data.

    Args:
        x (np.ndarray): Input 3D volume comprising of axial stacks of MRI images.

    Returns:
        np.ndarray: Mean-projected 2D projection in the coronal dimension.
    """
    return x.mean(axis=0).T


def project_sagittal(x: np.ndarray) -> np.ndarray:
    """Computes the mean 2D projection in the putative coronal dimension
    given axial input data.

    Args:
        x (np.ndarray): Input 3D volume comprising of axial stacks of MRI images.

    Returns:
        np.ndarray: Mean-projected 2D projection in the coronal dimension.
    """
    return x.mean(axis=1).T


def normalize(projection: np.ndarray) -> np.ndarray:
    """Normalize intensities according to the top-50 most intense
    voxels.

    Args:
        projection (np.ndarray): Input 2D projection.

    Returns:
        np.ndarray: Normalized 2D projection.
    """
    projection = 255 * projection / np.sort(projection)[:-50].mean()
    return projection.astype(np.uint16)


def center_pad(x: np.ndarray, width: int) -> np.ndarray:
    """Pad an image on the left and right with 0s to a specified
    target width.

    Args:
        x (np.ndarray): Input data.
        width (int): Desired width.

    Returns:
        np.ndarray: Padded data.
    """
    new_x = np.zeros((x.shape[0], width))
    offset = (width - x.shape[1]) // 2
    new_x[:, offset : width - offset] = x
    return new_x


def center_pad_stack(xs: np.ndarray) -> np.ndarray:
    """Center and pad input data and then stack the images
    with different widths.

    Args:
        xs (np.ndarray): Input data.

    Returns:
        np.ndarray: Center-padded and stacked data.
    """
    max_width = max(x.shape[1] for x in xs)
    return np.vstack([center_pad(x, max_width) for x in xs])


def build_z_slices(
    num_slices: List[int],
    z_pos: List[Tuple[float, float]],
) -> List[slice]:
    """
    Finds which images to take from each station.
    Takes a few noisy images from the top of each station
    and removes any remaining overlap between stations.

    num_slices: number of slices of each station
    z_pos: min z, max z of each station

    Example calculating the top and bottom z position of each station:
        z_pos = [50, 100], [90, 140]
        -> z_pos = [54, 100], [94, 140]  # remove images from the top
        -> overlap = [94, 100]
        -> z_pos = [54, 100], [100, 140]  # remove overlapping images
    """
    # remove a few slices from the top of each station
    top_remove_num = 4
    slice_starts = []
    for i in range(len(num_slices)):
        slice_starts.append(top_remove_num)
        top_remove_frac = top_remove_num / num_slices[i]
        lo, hi = z_pos[i]
        z_pos[i] = lo, hi - (hi - lo) * top_remove_frac

    # remove remaining overlaps from the bottom of each station
    slice_ends = []
    for i in range(0, len(num_slices) - 1):
        _, series_below_max_z = z_pos[i + 1]
        series_above_min_z, series_above_max_z = z_pos[i]
        overlap_size = series_below_max_z - series_above_min_z
        overlap_frac = overlap_size / (series_above_max_z - series_above_min_z)
        slice_ends.append(int((1 - overlap_frac) * num_slices[i]))
    slice_ends.append(
        num_slices[-1],
    )  # the last station gets nothing removed from the bottom
    # build slices
    return [slice(start, end) for start, end in zip(slice_starts, slice_ends)]


def build_projections(
    data: Dict[int, np.ndarray], meta_data: pd.DataFrame,
) -> Dict[str, np.ndarray]:
    """Build coronal and sagittal projections for each series type from all of the series.

    Args:
        data (Dict[int, np.ndarray]): Input data in the form {series number: series array}.
        meta_data (pd.DataFrame): Meta data from the Parquet files.

    Returns:
        Dict[str, np.ndarray]: Returns the dictionary {series type: projection}.
    """
    # Stations are differently scaled on the z-axis.
    station_z_scales = 3.0, 4.5, 4.5, 4.5, 3.5, 4.0
    station_z_scales = [scale / 3 for scale in station_z_scales]

    z_pos = meta_data.groupby("series_number")["image_position_z"].agg(["min", "max"])
    slices = build_z_slices(
        [data[i].shape[-1] for i in range(1, 25, 4)],
        [z_pos.loc[i] for i in range(1, 25, 4)],
    )

    # Keep track of where stations are connected by storing their intersection points in the
    # HDF5 dataset `horizontal_line_idx`.
    horizontal_lines = [
        (idx.stop - idx.start) * scale for idx, scale in zip(slices, station_z_scales)
    ]
    horizontal_lines = np.cumsum(horizontal_lines).astype(np.uint16)[:-1]
    projections = {"horizontal_line_idx": horizontal_lines}

    # Build coronal and sagittal projections.
    for type_idx, series_type_name in zip(range(4), ("in", "opp", "f", "w")):
        coronal_to_stack = []
        sagittal_to_stack = []
        for station_idx in range(1, 25, 4):  # neck, upper ab, lower ab, legs
            series_num = station_idx + type_idx
            station_slice = slices[station_idx // 4]
            scale = station_z_scales[station_idx // 4]
            coronal = project_coronal(data[series_num][..., station_slice])
            coronal = zoom(coronal, (scale, 1.0), order=1)  # account for z axis scaling
            coronal_to_stack.append(coronal)
            sagittal = project_sagittal(data[series_num][..., station_slice])
            sagittal = zoom(
                sagittal, (scale, 1.0), order=1,
            )  # Account for z axis scaling
            sagittal_to_stack.append(sagittal)

        projections[f"{series_type_name}_coronal"] = normalize(np.vstack(coronal_to_stack))
        projections[f"{series_type_name}_sagittal"] = normalize(
            center_pad_stack(sagittal_to_stack),
        )
    return projections


def build_projection_hd5(
        old_hd5_path: str,
        pq_base_path: str,
        output_folder: str,
):
    """Builds hd5 with 2d projections

    Args:
        old_hd5_path (str): Existing HDF5-file MRI slices.
        pq_base_path (str): Folder of existing meta data Parquet files.
        output_folder (str): Output path.
    """
    new_path = os.path.join(output_folder, os.path.basename(old_hd5_path))
    sample_id = os.path.splitext(os.path.basename(old_hd5_path))[0]
    with h5py.File(old_hd5_path, 'r') as old_hd5:
        for instance in old_hd5['instance']:
            data = {
                int(name): read_compressed(old_hd5[f'instance/{instance}/series/{name}'])
                for name in old_hd5[f'instance/{instance}/series']
            }
            meta_path = os.path.join(pq_base_path, f'{sample_id}_{instance}.pq')
            meta = ParquetFile(meta_path).to_pandas()
            projection = build_projections(data, meta)
            with h5py.File(new_path, 'a') as new_hd5:
                for name, im in projection.items():
                    compress_and_store(new_hd5, im, f'instance/{instance}/{name}')


def _build_projection_hd5s(hd5_files: List[str], pq_base_path: str, destination: str):
    """
    Applies build_projection_hd5 to a list of hd5 file paths and keeps track of errors.
    """
    errors = {}
    name = os.getpid()
    print(f'Starting process {name} with {len(hd5_files)} files')
    for i, path in enumerate(hd5_files):
        try:
            build_projection_hd5(path, pq_base_path, destination)
        except Exception as e:
            errors[path] = str(e)
        if len(hd5_files) % max(i // 10, 1) == 0:
            print(f'{name}: {(i + 1) / len(hd5_files):.2%} done')
    return errors


def multiprocess_project(
    hd5_files: List[str],
    pq_base_path: str,
    destination: str,
):
    """Builds hd5 with 2d projections

    Args:
        hd5_files (str): Existing HDF5-files containing MRI slices.
        pq_base_path (str): Folder of existing meta data Parquet files.
        destination (str): Output path.
    """
    os.makedirs(destination, exist_ok=True)
    split_files = np.array_split(hd5_files, cpu_count())
    print(f'Beginning coronal and sagittal projection of {len(hd5_files)} samples.')
    start = time.time()
    errors = {}
    with Pool(cpu_count()) as pool:
        results = [pool.apply_async(_build_projection_hd5s, (split, pq_base_path, destination)) for split in split_files]
        for result in results:
            errors.update(result.get())
    delta = time.time() - start
    print(f'Projections took {delta:.1f} seconds at {delta / len(hd5_files):.1f} s/file')
    with open(os.path.join(destination, 'errors.json'), 'w') as f:
        json.dump(errors, f)
    return errors
