import json
import logging
import numpy as np
import os
import pyro
import torch

from .cnv_model import SVDepthPyroModel, SVDepthData
from . import cnv_io


def run(args: dict, default_dtype: torch.dtype = torch.float32):
    base_path = os.path.join(args['model_dir'], args['model_name'])
    log_path = base_path + ".log.infer.txt"
    logging.basicConfig(filename=log_path,
                        filemode='w',
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    pyro.enable_validation(True)
    pyro.distributions.enable_validation(True)
    pyro.clear_param_store()

    np.random.seed(args['random_seed'])
    torch.random.manual_seed(args['random_seed'])
    pyro.set_rng_seed(args['random_seed'])

    params = load_model(base_path)
    if params is None:
        raise RuntimeError("Model at not found: {:s}".format(base_path))

    model = SVDepthPyroModel(k=params['k'], tensor_dtype=default_dtype, sample_depth_bin_size=params['sample_depth_bin_size'],
                             var_phi_bin=params['var_phi_bin'], var_phi_sample=params['var_phi_sample'],
                             alpha_non_ref=params['alpha_non_ref'], alpha_ref=params['alpha_ref'],
                             mu_eps=params['mu_eps'], device=args['device'], loss=params['loss'])
    load_param_store(base_path, device=args['device'])
    data = cnv_io.load_tensors(base_path=base_path, tensor_dtype=default_dtype, device=args['device'])
    logging.info(str(data.starts))

    posterior = model.run_predictive(data=data, n_samples=args['infer_predictive_samples'],
                                     n_iter=args['infer_predictive_iter'])

    discrete_posterior = model.run_discrete(data=data, n_states=model.k,
                                            log_freq=args['infer_discrete_log_freq'],
                                            n_samples=args['infer_discrete_samples'])
    posterior.update(discrete_posterior)

    # TODO get phi_sample
    output = get_output(data=data, stats=posterior, params=params)
    cnv_io.write_variant_output(output_path=args['output'], output_data=output)
    return output


def convert_type(dat):
    if isinstance(dat, np.ndarray):
        return dat.tolist()
    return dat


def get_output(data: SVDepthData, stats: dict, params: dict):
    n_bins = len(data.contigs)
    output_list = []
    for i in range(n_bins):
        output_list.append({
            'contig': data.contigs[i],
            'start': int(data.starts[i]),
            'bin_size': int(data.bin_size[i]),
            'freq_z': stats['z']['mean'][i, :],
            'ploidy': data.sample_ploidy[i, :].int().tolist(),
            'p_hw_loss': stats['p_hw']['mean'][i, 0],
            'p_hw_gain': stats['p_hw']['mean'][i, 2],
            'eps': stats['eps']['mean'][i] * params['mu_eps'],
            'phi_bin': stats['phi_bin']['mean'][i]
        })
    return output_list



def get_predictive_stats(samples: dict):
    return {key: {'mean': samples[key].astype(dtype='float').mean(axis=0).squeeze(),
                  'std': samples[key].astype(dtype='float').std(axis=0).squeeze()} for key in samples}


def get_discrete_stats(samples: dict):
    discrete_sites = ['m_pe', 'm_sr1', 'm_sr2'] #, 'm_rd']
    return {key: {'mean': samples[key].astype(dtype='float').mean(axis=0).squeeze()} for key in discrete_sites}


def calculate_state_frequencies(model: SVDepthPyroModel, discrete_samples: dict):
    z = discrete_samples['z']
    z_freq = np.zeros((z.shape[1], z.shape[2], model.k))
    for i in range(model.k):
        locs = z == i
        z_copy = z.copy()
        z_copy[locs] = 1
        z_copy[~locs] = 0
        z_freq[..., i] = z_copy.sum(axis=0)
    z_freq = z_freq / z_freq.sum(axis=2, dtype='float32', keepdims=True)
    m_pe_freq = discrete_samples['m_pe'].mean(axis=0)
    m_sr1_freq = discrete_samples['m_sr1'].mean(axis=0)
    m_sr2_freq = discrete_samples['m_sr2'].mean(axis=0)
    return {
        'z': z_freq,
        'm_pe': m_pe_freq,
        'm_sr1': m_sr1_freq,
        'm_sr2': m_sr2_freq
    }


def load_model(base_path: str):
    model_path = base_path + '.vars.json'
    if os.path.exists(model_path):
        with open(model_path, 'r') as f:
            return json.load(f)
    else:
        logging.warning("Could not locate model file {:s}".format(model_path))


def load_param_store(base_path: str, device: str = 'cpu'):
    pyro.clear_param_store()
    pyro.get_param_store().load(base_path + '.param_store.pyro', map_location=torch.device(device))
