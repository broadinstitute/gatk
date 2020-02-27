import json
import logging
import numpy as np
import os
import pyro
import torch

from . import constants
from .constants import SVTypes
from .model import SVGenotyperPyroModel
from . import io


def run(args):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    pyro.enable_validation(True)
    pyro.distributions.enable_validation(True)
    pyro.clear_param_store()

    np.random.seed(args.random_seed)
    torch.random.manual_seed(args.random_seed)
    pyro.set_rng_seed(args.random_seed)

    output_by_type = {}
    global_stats_by_type = {}
    for svtype in [SVTypes.DEL, SVTypes.DUP, SVTypes.INS, SVTypes.INV]:
        base_path = os.path.join(args.model_dir, args.model_name + "." + svtype.name)
        params = load_model(base_path)
        if params is None:
            logging.info("Skipping SV type {:s}".format(svtype.name))
        else:
            model = SVGenotyperPyroModel(svtype=svtype, k=params['k'], mu_eps_pe=params['mu_eps_pe'], mu_eps_sr1=params['mu_eps_sr1'],
                                         mu_eps_sr2=params['mu_eps_sr2'], mu_lambda_pe=params['mu_lambda_pe'], mu_lambda_sr1=params['mu_lambda_sr1'],
                                         mu_lambda_sr2=params['mu_lambda_sr2'], var_phi_pe=params['var_phi_pe'], var_phi_sr1=params['var_phi_sr1'],
                                         var_phi_sr2=params['var_phi_sr2'], mu_eta_q=params['mu_eta_q'], mu_eta_r=params['mu_eta_r'],
                                         device=args.device, loss=params['loss'])
            load_param_store(base_path, device=args.device)
            vids_list = io.load_list(base_path + ".vids.list")
            sample_ids_list = io.load_list(base_path + ".sample_ids.list")
            data = io.load_tensors(directory=args.model_dir, model_name=args.model_name, svtype=svtype, device=args.device)

            predictive_samples = model.infer_predictive(data=data, n_samples=args.infer_predictive_samples)
            discrete_samples = model.infer_discrete(data=data, svtype=svtype, log_freq=args.infer_discrete_log_freq, n_samples=args.infer_discrete_samples)
            freq = calculate_state_frequencies(model=model, discrete_samples=discrete_samples)
            genotypes = get_genotypes(freq_z=freq['z'])
            stats = get_predictive_stats(samples=predictive_samples)
            stats.update(get_discrete_stats(samples=discrete_samples))
            output_by_type[svtype] = get_output(vids_list=vids_list, genotypes=genotypes, stats=stats, params=params)
            global_stats_by_type[svtype] = get_global_stats(stats)
    output = {}
    for svtype in output_by_type:
        output.update(output_by_type[svtype])
    return output, global_stats_by_type


def get_output(vids_list: list, genotypes: dict, stats: dict, params: dict):
    n_variants = len(vids_list)
    output_dict = {}
    for i in range(n_variants):
        vid = vids_list[i]
        if 'eps_pe' in stats:
            eps_pe = stats['eps_pe']['mean'][i] * params['mu_eps_pe']
        else:
            eps_pe = 0
        if 'phi_pe' in stats:
            phi_pe = stats['phi_pe']['mean'][i]
        else:
            phi_pe = 0
        if 'eta_r' in stats:
            eta_r = stats['eta_r']['mean'][i] * params['mu_eta_r']
        else:
            eta_r = 0
        output_dict[vid] = {
            'gt': genotypes['gt'][i, :],
            'gt_p': genotypes['gt_p'][i, :],
            'gt_lod': genotypes['gt_lod'][i, :],
            'p_m_pe': stats['m_pe']['mean'][i],
            'p_m_sr1': stats['m_sr1']['mean'][i],
            'p_m_sr2': stats['m_sr2']['mean'][i],
            'p_m_rd': stats['m_rd']['mean'][i],
            'eps_pe': eps_pe,
            'eps_sr1': stats['eps_sr1']['mean'][i] * params['mu_eps_sr1'],
            'eps_sr2': stats['eps_sr2']['mean'][i] * params['mu_eps_sr2'],
            'phi_pe': phi_pe,
            'phi_sr1': stats['phi_sr1']['mean'][i],
            'phi_sr2': stats['phi_sr2']['mean'][i],
            'eta_q': stats['eta_q']['mean'][i] * params['mu_eta_q'],
            'eta_r': eta_r
        }
    return output_dict


def get_global_stats(stats: dict):
    global_sites = ['pi_sr1', 'pi_sr2', 'pi_pe', 'pi_rd', 'lambda_pe', 'lambda_sr1', 'lambda_sr2']
    global_stats = {}
    for site in global_sites:
        if site in stats:
            global_stats[site + '_mean'] = stats[site]['mean']
            global_stats[site + '_std'] = stats[site]['std']
    return global_stats


def get_predictive_stats(samples: dict):
    return {key: {'mean': samples[key].astype(dtype='float').mean(axis=0).squeeze(),
                  'std': samples[key].astype(dtype='float').std(axis=0).squeeze()} for key in samples}


def get_discrete_stats(samples: dict):
    discrete_sites = ['m_pe', 'm_sr1', 'm_sr2', 'm_rd']
    return {key: {'mean': samples[key].astype(dtype='float').mean(axis=0).squeeze()} for key in discrete_sites}


def get_genotypes(freq_z: dict):
    gt = freq_z.argmax(axis=2)
    freq_sorted = np.sort(freq_z, axis=2)
    gt_p = freq_sorted[..., -1]
    gt_lod = np.log(freq_sorted[..., -1]) - np.log(freq_sorted[..., -2])
    gt_lod[np.isinf(gt_lod)] = constants.MAX_GT_LOD
    gt_lod[gt_lod > constants.MAX_GT_LOD] = constants.MAX_GT_LOD
    return {
        'gt': gt,
        'gt_p': gt_p,
        'gt_lod': gt_lod
    }


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
