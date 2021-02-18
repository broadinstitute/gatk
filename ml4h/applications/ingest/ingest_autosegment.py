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
import os
import h5py
import numpy as np
import pandas as pd
from scipy.ndimage import zoom
from fastparquet import ParquetFile
import cv2
import blosc
from ingest_mri import compress_and_store
from two_d_projection import build_z_slices


def uncompress(t: np.ndarray, stored_dtype=np.uint16):
    return np.frombuffer(blosc.decompress(t[()]), dtype=stored_dtype).reshape(
        t.attrs["shape"]
    )


def pad_center(img: np.ndarray, shape):
    border_v = 0
    border_h = 0
    if (shape[0] / shape[1]) >= (img.shape[0] / img.shape[1]):
        border_v = int((((shape[0] / shape[1]) * img.shape[1]) - img.shape[0]) / 2)
    else:
        border_h = int((((shape[1] / shape[0]) * img.shape[0]) - img.shape[1]) / 2)
    img = cv2.copyMakeBorder(
        img, border_v, border_v, border_h, border_h, cv2.BORDER_CONSTANT, 0
    )
    img = cv2.resize(img, (shape[1], shape[0]))
    return img


def autosegment_axial_slice(img: np.ndarray):
    """Given an input axial slice we will try to automatically segment
    its contours using standard image processing techniques.

    Args:
        img (np.ndarray): Input axial slice.

    Returns:
        closing: 2D segmentation of dimensions equal to the input axial slice
        is_legs: Boolean flag for whether or not the input axial slice likely contains legs
        contour_length: Arc-length (circumference) of the 2D segmentation in pixels
    """
    img = (img / img.max() * 255).astype(np.uint8)
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    img = clahe.apply(img)
    _, thresh = cv2.threshold(img, 0, 255, cv2.THRESH_OTSU)
    _, markers = cv2.connectedComponents(thresh)
    marker_area = [np.sum(markers == m) for m in range(np.max(markers) + 1) if m != 0]
    marker_area_rank = np.argsort(marker_area)[::-1]
    top2 = np.array(marker_area)[marker_area_rank[:2]]

    countour_length = 0
    is_legs = False
    if ((top2 / top2.max()).min() >= 0.25) and (len(top2) == 2):
        is_legs = True
        thresh = (
            (markers == (marker_area_rank + 1)[0])
            | (markers == (marker_area_rank + 1)[1])
        ).astype(np.uint8)
    else:
        largest_component = np.argmax(marker_area) + 1
        thresh = (markers == largest_component).astype(np.uint8)

    # Fill holes
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    closing = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, k, iterations=1)
    contour, _ = cv2.findContours(closing, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    for cnt in contour:
        cv2.drawContours(closing, [cnt], 0, 255, -1)

    # Recapture surface area
    contour, _ = cv2.findContours(closing, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    for cnt in contour:
        countour_length += cv2.arcLength(cnt, True)

    return closing, is_legs, countour_length


def autosegment_dixon(meta: str, file: str, destination: str, instance: int = 2):
    """Given an input HDF5 file of data and Parquet file of meta data, we will
    try to automatically segment its contours using standard image processing techniques.

    This subroutine includes several constants and hardcoded hyperparameters
    specific to the UK Biobank whole-body Dixon MRI data. No guarantees are
    provided for other related data.

    Args:
        meta (str): Path to Parquet file with DICOM meta data on disk
        file (str): Path to HDF5 file with data on disk
        destination (str): Output path on disk
        instance (int, optional): UK Biobank instance number. Defaults to 2.

    Returns:
        Returns True
    """
    meta = ParquetFile(meta).to_pandas()

    with h5py.File(file, "r") as f:
        data = {
            int(name): uncompress(f[f"instance/{instance}/series/{name}"])
            for name in f[f"instance/{instance}/series"]
        }

    z_pos = meta.groupby("series_number")["image_position_z"].agg(["min", "max"])
    slices = build_z_slices(
        [data[i].shape[-1] for i in range(1, 25, 4)],
        [z_pos.loc[i] for i in range(1, 25, 4)],
    )

    # Sample data to 3mm resolution.
    total = []
    scaling_factors = [3.0, 4.5, 4.5, 4.5, 3.5, 4.0]
    for i, s in zip(range(1, 25, 4), scaling_factors):
        d = data[i][..., slices[i // 4].start : slices[i // 4].stop]
        d = zoom(d, (1.0, 1.0, s / 3.0))
        d = pad_center(d, (174, 224, d.shape[-1]))
        total.append(d)

    dat = np.concatenate(total, axis=-1)

    # X,Y,Z dimension resolution in millimeters (mm).
    xmm = meta["col_pixel_spacing_mm"].iloc[0]  # 2.232142925262451 mm
    ymm = meta["row_pixel_spacing_mm"].iloc[0]  # 2.232142925262451 mm
    zmm = 3.0  # mm

    # Store axial data in arrays.
    stack = []
    leg_flag = []
    cubic_mm = []
    surface_area = []

    # Iterate over axial slices.
    for i in range(dat.shape[-1]):
        closing, is_legs, clen = autosegment_axial_slice(dat[..., i])
        stack.append(closing)
        leg_flag.append(is_legs)
        cubic_mm.append((closing == 255).sum() * (xmm * ymm * zmm))
        surface_area.append(clen)

    # Stack axial contour slices in 3D volume and cast to
    # 2-byte for immediate compatibility with existing ingest
    # code.
    stack = np.array(stack).astype(np.uint16)

    # Dataframe of volumes
    p = pd.DataFrame(
        {
            "x": range(len(cubic_mm)),  # Number of axial slices
            "volume": cubic_mm,  # Volume in mm3
            "surface_area": surface_area,  # Contour length in pixels
            "is_leg": leg_flag,  # Boolean flag for whether we believe this axial slice is a leg
        }
    )
    p["volume"] /= 1e6  # Convert to L
    p["surface_area"] *= xmm * zmm  # pixel length * mm/pixel * depth of stack in mm
    p["surface_area"] /= 1e6  # Convert to mm2
    p["x_rev"] = p["x"][::-1].values  # Also store reverse range for plotting reasons
    p = p.loc[p.index.values[::-1]]  # Store in reverse order

    # Store computed results
    output_name = str(meta.ukbid.iloc[0])
    p.to_parquet(os.path.join(destination, f"{output_name + '.pq'}"))

    # Store 3D surface and axial stack
    with h5py.File(os.path.join(destination, f"{output_name + '.h5'}"), "a") as f:
        # Save the 3D surface stack.
        compress_and_store(f, stack, f"/instance/{instance}/series/surface")
        # Save the 3D axial stack used for the surface.
        compress_and_store(f, dat, f"/instance/{instance}/series/axial_stack")

    return True


def _build_autosegment_dixon_hd5s(df: pd.DataFrame, destination: str):
    errors = {}
    name = os.getpid()
    print(f"Starting process {name} with {len(df)} files")

    for i in range(len(df)):
        try:
            autosegment_dixon(df.meta.iloc[i], df.file.iloc[i], destination)
        except Exception as e:
            errors[df.file.iloc[i]] = str(e)
        if len(df) % max(i // 10, 1) == 0:
            print(f"{name}: {(i + 1) / len(df):.2%} done")

    return errors


def multiprocess_autosegment_dixon(
    df: pd.DataFrame,
    destination: str,
):
    """Ingest autosegmentations of Dixon volumes. Takes a DataFrame with two columns
    as input:
        `file`: String path to HDF5 file on disk
        `meta`: String path to a Parquet file with DICOM meta data on disk.

    Args:
        df (pd.DataFrame): Input DataFrame with file/meta paths.
        destination (str): Output destination on disk

    Returns:
        dict: A dictionary of caught exceptions
    """
    os.makedirs(destination, exist_ok=True)
    split_files = np.array_split(df, cpu_count())
    print(f"Beginning autosegmentation of {len(df)} samples.")
    start = time.time()
    errors = {}

    with Pool(cpu_count()) as pool:
        results = [
            pool.apply_async(_build_autosegment_dixon_hd5s, (split, destination))
            for split in split_files
        ]
        for result in results:
            errors.update(result.get())

    delta = time.time() - start
    print(f"Autosegmentation took {delta:.1f} seconds at {delta / len(df):.1f} s/file")

    with open(os.path.join(destination, "errors.json"), "w") as f:
        json.dump(errors, f)

    return errors
