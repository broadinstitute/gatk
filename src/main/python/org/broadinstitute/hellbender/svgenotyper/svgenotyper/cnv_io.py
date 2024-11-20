import numpy as np
import pandas as pd
import torch

from . import constants
from .cnv_model import SVDepthData


def load_batch(variants_file_path: str,
               device: str,
               tensor_dtype: torch.dtype):
    contig_list = []
    start_list = []
    bin_size_list = []
    counts_list = []
    ploidy_list = []
    with open(variants_file_path, 'r') as f:
        for line in f:
            fifo_data = line.split('\t')
            contig_list.append(fifo_data[0])
            start_list.append(int(fifo_data[1]))
            bin_size_list.append(int(fifo_data[2]))
            counts_list.append([int(x) for x in fifo_data[3].split(';')])
            ploidy_list.append([int(x) for x in fifo_data[4].split(';')])
    contig_np = np.asarray(contig_list)
    start_t = torch.tensor(start_list, device=device, dtype=tensor_dtype)
    bin_size_t = torch.tensor(bin_size_list, device=device, dtype=tensor_dtype)
    counts_t = torch.tensor(counts_list, device=device, dtype=tensor_dtype)
    ploidy_t = torch.tensor(ploidy_list, device=device, dtype=tensor_dtype)
    return contig_np, start_t, bin_size_t, counts_t, ploidy_t


def load_data(variants_file_path: str,
              mean_coverage_path: str,
              samples_path: str,
              tensor_dtype: torch.dtype,
              device: str = 'cpu'):
    # TODO: cross-check depth table samples with samples list
    sample_count_df = pd.read_csv(mean_coverage_path, sep='\t', header=None, index_col=0)
    sample_depth_t = torch.from_numpy(sample_count_df.values).to(device=device, dtype=tensor_dtype).squeeze(-1) / torch.tensor(constants.DEPTH_PLOIDY).to(device=device, dtype=tensor_dtype)
    samples_np = np.loadtxt(samples_path, dtype=str)
    contig_np, start_t, bin_size_t, counts_t, ploidy_t = load_batch(variants_file_path=variants_file_path, device=device, tensor_dtype=tensor_dtype)
    if contig_np.shape[0] == 0:
        return None
    return SVDepthData(samples=samples_np, sample_per_base_depth=sample_depth_t, sample_ploidy=ploidy_t, contigs=contig_np, starts=start_t, bin_size=bin_size_t, counts=counts_t)


def write_variant_output(output_path: str, output_data: dict):
    param_keys = ["eps", "phi_bin", "p_hw_loss", "p_hw_gain"]
    header = ["contig", "start", "bin_size", "freq_z"] + param_keys
    with open(output_path, 'w') as f, open("temp", 'w') as fd:
        line = "#" + "\t".join(header)
        f.write(line + "\n")
        for d in output_data:
            z_freq = pretty_print_2d_array(d['freq_z'])
            line = "\t".join([d['contig'], str(d['start']), str(d['bin_size']), z_freq] + [str(d[x]) for x in param_keys])
            f.write(line + "\n")
            fd.write(line + "\n")


def pretty_print_2d_array(arr):
    return ";".join(",".join(str(y) for y in x) for x in arr.tolist())


def save_tensors(data: SVDepthData, base_path: str):
    data_vars = vars(data)
    for var in data_vars:
        if data_vars[var] is not None:
            torch.save(data_vars[var], base_path + "." + var + ".pt")


def save_list(data: list, path: str):
    with open(path, 'w') as f:
        f.writelines([x + '\n' for x in data])


def load_tensors(base_path: str, tensor_dtype: torch.dtype, device: str = 'cpu'):
    sample_per_base_depth = torch.load(base_path + ".sample_per_base_depth.pt", map_location=device).to(dtype=tensor_dtype)
    sample_ploidy = torch.load(base_path + ".sample_ploidy.pt", map_location=device).to(dtype=tensor_dtype)
    starts = torch.load(base_path + ".starts.pt", map_location=device).to(dtype=tensor_dtype)
    bin_size = torch.load(base_path + ".bin_size.pt", map_location=device).to(dtype=tensor_dtype)
    counts = torch.load(base_path + ".counts.pt", map_location=device).to(dtype=tensor_dtype)
    samples = np.loadtxt(base_path + ".sample_ids.list", dtype=str)
    contigs = np.loadtxt(base_path + ".contigs.list", dtype=str)
    return SVDepthData(samples=samples, sample_per_base_depth=sample_per_base_depth, sample_ploidy=sample_ploidy,
                       contigs=contigs, starts=starts, bin_size=bin_size, counts=counts)


def load_list(path: str):
    with open(path, 'r') as f:
        lines = f.readlines()
        if lines is None or len(lines) == 0:
            raise ValueError("Tried to read empty list: {:s}".format(path))
        return [x.strip() for x in lines]
