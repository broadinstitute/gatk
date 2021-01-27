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
import ml4h
import tensorflow as tf
import pandas as pd
from tensorflow.keras.models import load_model
from tensorflow_addons.optimizers import RectifiedAdam
import numpy as np
import h5py
import os
import glob
from ml4h.metrics import get_metric_dict
import ml4h.tensormap.ukb.mri
from ml4h.tensormap.general import pad_or_crop_array_to_shape
from ml4h.normalizer import ZeroMeanStd1
import zstandard
import blosc
import pyarrow
from io import BytesIO
import sys
import timeit


def prepare_model(model_file: str, tensormap: ml4h.tensormap) -> tf.keras.Model:
    """Load a pre-trained ML4H model from disk and all appropriate
    custom objects.

    Args:
        model_file (str): Path to TF/Keras model.
        tensormap (ml4h.tensormap): TensorMap used to train/evaluate the model.

    Returns:
        tf.keras.Model: Returns a Keras model.
    """
    # Get the config using the TensorMap used to train the model in the first place.
    # This allows us to reconstruct models that are not saved using configs.
    objects = get_metric_dict([tensormap])
    # Silly work-around to reset function pointer to our local instance.
    objects["RectifiedAdam"] = RectifiedAdam
    # Load the model
    model = load_model(model_file, custom_objects=objects)
    model.summary()
    return model


def split_files_for_parallel_computing(
    files, partition_number: int, total_partitions: int,
):
    """Given a list of files, return the Nth partition of files. This is
    function is used when distribution the workload across multiple machines
    in an embarrassingly parallel fashion.

    Args:
        files ([str]): Input list of file paths.
        partition_number (int): The target batch/partition.
        total_partitions (int): Total number of desired subpartitions.

    Returns:
        [str]: Returns the subpartition of files that corresponds to this partition number
    """
    step = len(files) // total_partitions
    if partition_number != step:
        files = files[(step * (partition_number)) : (step * (partition_number + 1))]
    else:
        files = files[(step * (partition_number)) :]
    return files


def prepare_local_tensor(i: int, tensor, names):
    """Special tensorization callback function for the model
    `sax_slices_jamesp_4b_hyperopted_dropout_pap_dupe`.
    """
    tensor_local = np.zeros((50, 224, 224, 4))
    # Uncomment for debugging:
    # print(max(i - 2, 0), max(i - 1, 0), i, min(i + 1, len(names) - 1))
    if i - 2 >= 0:
        tensor_local[:, :, :, 0] = tensor[(i - 2), ...]
    else:
        if -1 >= 0:
            tensor_local[:, :, :, 0] = tensor[(i - 1), ...]
        else:
            tensor_local[:, :, :, 0] = tensor[(i), ...]
    if i - 1 >= 0:
        tensor_local[:, :, :, 1] = tensor[(i - 1), ...]
    else:
        tensor_local[:, :, :, 1] = tensor[0, ...]
    tensor_local[:, :, :, 2] = tensor[(i), ...]
    if i + 1 < len(names):
        tensor_local[:, :, :, 3] = tensor[(i + 1), ...]
    else:
        tensor_local[:, :, :, 3] = tensor[i, ...]
    return tensor_local


def jpp_infer_short_axis(files, model: tf.keras.Model, output_path: str):
    """Inference loop that takes a list of prepared HDF5 files,
    a pre-trained model, and an output path, and computes the
    channel-wise argmax of the inference result.

    Args:
        files ([type]): [description]
        model (tf.keras.Model): [description]
        output_path (str): [description]
    """
    tot = 0
    for s in files:
        print(f"{s}. Progress: {tot}/{len(files)}")
        try:
            x = h5py.File(s, "r")
        except Exception as e:
            print(f"Failed to open {s}")
            continue

        try:
            names = pd.DataFrame({"path": list(x["/ukb_cardiac_mri/"])})
        except Exception as e:
            print(f"Malformed hdf5 {s}")
            continue

        names = names[names.loc[:, "path"].str.contains("cine_segmented_sax_b")]
        names["slices"] = [
            int(s.split("_")[-1].replace("b", "")) for s in names["path"].values
        ]
        names = names.sort_values(["slices"])
        names = names[~names.path.str.contains("james")]

        if len(names) == 0:
            print(f"Failed to find data for {s}")
            continue

        instances = list(x["/ukb_cardiac_mri/"][names.path.iloc[0]])

        for instance in instances:
            print(f"instance: {instance}")
            names_valid = names[
                [instance in list(x["/ukb_cardiac_mri/"][p]) for p in names.path]
            ]
            # Grab data
            tensor = np.zeros((len(names_valid), 50, 224, 224), dtype=np.float32)
            for f, k in zip(names_valid.path, range(len(names_valid))):
                tensor[k, ...] = ZeroMeanStd1().normalize(
                    pad_or_crop_array_to_shape(
                        (50, 224, 224),
                        np.moveaxis(
                            x["/ukb_cardiac_mri/"][f][instance]["instance_0"][()], 2, 0,
                        ),
                    ),
                )

            start_predict = timeit.default_timer()  # Debug timer
            test = np.zeros((len(names_valid), 50, 224, 224, 17), dtype=np.float32)
            # These are no longer used. Left for legacy reasons.
            # argmax = np.zeros((len(names_valid), 50, 224, 224),     dtype=np.float32)

            for i in range(len(names_valid)):
                tensor_local = prepare_local_tensor(i, tensor, names_valid)
                for j in range(50):
                    pred = model.predict(tensor_local[j : (j + 1), ...])
                    test[i, j, ...] = tf.squeeze(pred)
                    # These are no longer used. Left for legacy reasons.
                    # argmax[i,j, ...] = tf.argmax(pred, axis=-1)

            # Flow condition: argmax or prob
            filename = os.path.splitext(os.path.split(s)[-1])[0]
            with h5py.File(
                os.path.join(output_path, f"{filename}_inference__argmax.h5"), "a",
            ) as ff:
                ff_in = ff.create_group(f"instance_{str(instance)}")
                ff_in.create_dataset(
                    "argmax",
                    data=np.void(
                        blosc.compress(
                            tf.argmax(test, axis=-1).numpy().astype(np.uint8).tobytes(),
                            typesize=2,
                            cname="zstd",
                            clevel=9,
                        ),
                    ),
                )
                ff_in.attrs["shape"] = test.shape
                # Hard-core approach to store Parquet as an in-memory view
                buffer = BytesIO()
                names_valid.to_parquet(buffer, engine="pyarrow", compression="zstd")
                ff_in.create_dataset("slices_pq", data=np.void(buffer.getvalue()))
                # Getting data back:
                # pd.read_parquet(BytesIO(buffer), engine='pyarrow')

            stop_predict = timeit.default_timer()  # Debug timer
            print("Predict time: ", stop_predict - start_predict)  # Debug message
            tot += 1
