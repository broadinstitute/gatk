import json
import logging
import numpy as np
import os
import pyro
import torch

from .constants import SVTypes
from .model import SVGenotyperPyroModel
from . import io


def run(args: dict, default_dtype: torch.dtype = torch.float32):
    base_path = os.path.join(args['model_dir'], args['model_name'])
    log_path = base_path + ".log.genotype.txt"
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
    svtype = SVTypes[args['svtype']]

    params = load_model(base_path)
    if params is None:
        raise RuntimeError("Model at not found: {:s}".format(base_path))

    model = SVGenotyperPyroModel(svtype=svtype, k=params['k'], tensor_dtype=default_dtype, mu_eps_pe=params['mu_eps_pe'], mu_eps_sr1=params['mu_eps_sr1'],
                                 mu_eps_sr2=params['mu_eps_sr2'], var_phi_pe=params['var_phi_pe'], var_phi_sr1=params['var_phi_sr1'],
                                 var_phi_sr2=params['var_phi_sr2'], read_length=params['read_length'],
                                 device=args['device'], loss=params['loss'])
    load_param_store(base_path, device=args['device'])
    vids_list = io.load_list(base_path + ".vids.list")
    data = io.load_tensors(base_path=base_path, svtype=svtype, tensor_dtype=default_dtype, device=args['device'])

    posterior = model.get_quantiles(data=data)
    discrete_posterior = model.run_discrete(data=data, svtype=svtype, n_states=model.k,
                                            log_freq=args['genotype_discrete_log_freq'],
                                            n_samples=args['genotype_discrete_samples'])
    posterior.update(discrete_posterior)
    output = get_output(vids_list=vids_list, stats=posterior, params=params)
    io.write_variant_output(output_path=args['output'], output_data=output)
    return output


def convert_type(dat):
    if isinstance(dat, np.ndarray):
        return dat.tolist()
    return dat


def get_output(vids_list: list, stats: dict, params: dict):
    n_variants = len(vids_list)
    output_list = []
    for i in range(n_variants):
        if 'eps_pe' in stats:
            eps_pe_median = stats['eps_pe']['median'][i] * params['mu_eps_pe']
            eps_pe_iqr = stats['eps_pe']['iqr'][i] * params['mu_eps_pe']
        else:
            eps_pe_median = 0
            eps_pe_iqr = 0
        if 'phi_pe' in stats:
            phi_pe_median = stats['phi_pe']['median'][i]
            phi_pe_iqr = stats['phi_pe']['iqr'][i]
        else:
            phi_pe_median = 0
            phi_pe_iqr = 0
        if 'm_pe' in stats:
            p_m_pe_mean = stats['m_pe']['mean'][i]
        else:
            p_m_pe_mean = 0
        if 'q_hw' in stats:
            q_hw_median = stats['q_hw']['median'][i]
            q_hw_iqr = stats['q_hw']['iqr'][i]
        else:
            q_hw_median = 0
            q_hw_iqr = 0
        if 'r_hw' in stats:
            r_hw_median = stats['r_hw']['median'][i]
            r_hw_iqr = stats['r_hw']['iqr'][i]
        else:
            r_hw_median = 0
            r_hw_iqr = 0
        output_list.append({
            'vid': vids_list[i],
            'freq_z': stats['z']['mean'][i, :],
            'p_m_pe': p_m_pe_mean,
            'p_m_sr1': stats['m_sr1']['mean'][i],
            'p_m_sr2': stats['m_sr2']['mean'][i],
            'eps_pe_median': eps_pe_median,
            'eps_sr1_median': stats['eps_sr1']['median'][i] * params['mu_eps_sr1'],
            'eps_sr2_median': stats['eps_sr2']['median'][i] * params['mu_eps_sr2'],
            'phi_pe_median': phi_pe_median,
            'phi_sr1_median': stats['phi_sr1']['median'][i],
            'phi_sr2_median': stats['phi_sr2']['median'][i],
            'q_hw_median': q_hw_median,
            'r_hw_median': r_hw_median,
            'eps_pe_iqr': eps_pe_iqr,
            'eps_sr1_iqr': stats['eps_sr1']['iqr'][i] * params['mu_eps_sr1'],
            'eps_sr2_iqr': stats['eps_sr2']['iqr'][i] * params['mu_eps_sr2'],
            'phi_pe_iqr': phi_pe_iqr,
            'phi_sr1_iqr': stats['phi_sr1']['iqr'][i],
            'phi_sr2_iqr': stats['phi_sr2']['iqr'][i],
            'q_hw_iqr': q_hw_iqr,
            'r_hw_iqr': r_hw_iqr
        })
    return output_list


def get_predictive_stats(samples: dict):
    return {key: {'mean': samples[key].astype(dtype='float').mean(axis=0).squeeze(),
                  'std': samples[key].astype(dtype='float').std(axis=0).squeeze()} for key in samples}


def get_discrete_stats(samples: dict):
    discrete_sites = ['m_pe', 'm_sr1', 'm_sr2'] #, 'm_rd']
    return {key: {'mean': samples[key].astype(dtype='float').mean(axis=0).squeeze()} for key in discrete_sites}


def calculate_state_frequencies(model: SVGenotyperPyroModel, discrete_samples: dict):
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
