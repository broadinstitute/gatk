# Reproduction scripts: inference RV

This folder contains scripts to reproduce the models used in:

**Genetic Analysis of Right Heart Structure and Function in 45,000 People**. James P. Pirruccello*, Paolo Di Achille*, Victor Nauffal*, Mahan Nekoui, Samuel N. Friedman, Marcus D. R. Klarqvist, Mark D. Chaffin, Shaan Khurshid, Carolina Roselli, Puneet Batra, Kenney Ng, Steven A. Lubitz, Jennifer E. Ho, Mark E. Lindsay, Anthony Philippakis, Patrick T. Ellinor. [To appear]

## Example

Given a pre-trained semantic segmentation model `sax_slices_jamesp_4b_hyperopted_dropout_pap_dupe.h5` and the `ml4h.tensormap` that was used to generate data we can proceed to make inference on new data.

```py
import infer_on_sax # Local file

# Pre-trained model
model = prepare_model("/tf/sax_slices_jamesp_4b_hyperopted_dropout_pap_dupe.h5", ml4h.tensormap.ukb.mri.cine_segmented_sax_slice_jamesp)
# Enumerate the target files of interest.
files = glob.glob('/mnt/disks/annotated-cardiac-tensors-44k/2020-09-21/*.hd5')
# Partition the files into buckets and retrieve the files corresponding to that bucket.
# For example, embarassingly parallel computation across 50 GCP VMs with NVidia P4 GPUs
# using the provided shell script.
files = split_files_for_parallel_computing(files, partition_number=0, total_partitions=50)
jpp_infer_short_axis(files, model, output_path='/tf/')
```

A provided shell script `infer_hdf5_to_local.sh` streamline the procedure of spawning multiple GCP VMs with attached disks and GPUs for inference. Make sure you modify this file for executing the appropriate commands on the VMs.
