import numpy as np

import pyro
from pyro.infer import JitTraceEnum_ELBO, TraceEnum_ELBO
from pyro.infer.svi import SVI
from pyro.optim import PyroOptim

import torch
from torch.distributions import constraints

from . import constants
from .constants import SVTypes
from .model import SVGenotyperData, SVGenotyperPyroModel
from .io import load_data

import logging
from typing import List


def train_epoch(svi: SVI, data: SVGenotyperData, epoch: int, scheduler: pyro.optim.LambdaLR) -> float:
    loss = svi.step(data_pe=data.pe_t, data_sr1=data.sr1_t, data_sr2=data.sr2_t, depth_t=data.depth_t,
                    rd_gt_prob_t=data.rd_gt_prob_t)
    scheduler.step(epoch)
    return loss


def create_scheduler(lr_min: float, lr_init: float, lr_decay: float, beta1: float, beta2: float):
    return pyro.optim.LambdaLR({
        'optimizer': torch.optim.Adam,
        'optim_args': {'lr': 1., 'betas': (beta1, beta2)},
        'lr_lambda': lambda k: lr_min + (lr_init - lr_min) * np.exp(-k / lr_decay)
    })


def initialize_guide(svi: SVI, model: SVGenotyperPyroModel, data: SVGenotyperData) -> float:
    n_latent_dim = model.get_latent_dim(data.pe_t.shape[0])
    loc_init = torch.zeros(n_latent_dim, device=model.device)
    pyro.param("auto_loc", loc_init)
    scale_init = torch.ones(n_latent_dim, device=model.device)
    pyro.param("auto_scale", scale_init, constraint=constraints.positive)
    return svi.loss(model.model, model.guide, data_pe=data.pe_t, data_sr1=data.sr1_t, data_sr2=data.sr2_t,
                    depth_t=data.depth_t, rd_gt_prob_t=data.rd_gt_prob_t)


def run_training(model: SVGenotyperPyroModel,
                 data: SVGenotyperData,
                 args) -> List[float]:
    logging.info("Initializing model...")
    pyro.clear_param_store()
    scheduler = create_scheduler(lr_min=args.lr_min, lr_init=args.lr_init, lr_decay=args.lr_decay,
                                 beta1=args.adam_beta1, beta2=args.adam_beta2)
    if args.jit:
        loss = JitTraceEnum_ELBO()
    else:
        loss = TraceEnum_ELBO()
    svi = SVI(model.model, model.guide, optim=scheduler, loss=loss)
    initial_loss = initialize_guide(svi=svi, model=model, data=data)
    train_elbo = []

    logging.info("Initial loss: {:.4f}".format(initial_loss))
    logging.info("Running model training...")
    # Run training loop.  Use try to allow for keyboard interrupt.
    try:
        for epoch in range(args.max_iter):
            # Train, and keep track of training loss.
            total_epoch_loss = train_epoch(svi=svi, data=data, epoch=epoch, scheduler=scheduler)
            train_elbo.append(-total_epoch_loss)
            model.loss['train']['epoch'].append(epoch)
            model.loss['train']['elbo'].append(-total_epoch_loss)
            if (epoch + 1) % args.iter_log_freq == 0:
                logging.info("[epoch %04d]  training loss: %.4f" % (epoch + 1, total_epoch_loss))
        logging.info("Training procedure complete after {:d} epochs, loss: {:.4f}".format(args.max_iter, total_epoch_loss))

    # Exception allows program to continue after ending training prematurely.
    except KeyboardInterrupt:
        logging.info("Training procedure stopped by keyboard interrupt.")

    return train_elbo


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
    for svtype in [SVTypes.DEL, SVTypes.DUP, SVTypes.INS, SVTypes.INV]:
        model = SVGenotyperPyroModel(svtype, args.num_states, args.eps_pe, args.eps_sr1, args.eps_sr2, args.device)
        vids_np, data = load_data(vcf_path=args.vcf,
                                  mean_coverage_path=args.coverage_file,
                                  svtype=svtype,
                                  num_states=args.num_states,
                                  depth_dilution_factor=args.depth_dilution_factor,
                                  device=args.device)
        if data is None:
            logging.info("No records of type {:s} found.".format(str(svtype.name)))
        else:
            logging.info("Training {:s} with {:d} variants and {:d} samples...".format(str(svtype.name), data.pe_t.shape[0], data.pe_t.shape[1]))
            run_training(model=model, data=data, args=args)
            predictive_samples = model.infer_predictive(data=data, log_freq=args.infer_predictive_log_freq, n_samples=args.infer_predictive_samples)
            discrete_samples = model.infer_discrete(data=data, svtype=svtype, log_freq=args.infer_discrete_log_freq, n_samples=args.infer_discrete_samples)
            freq = calculate_state_frequencies(model=model, discrete_samples=discrete_samples)
            genotypes = get_genotypes(freq_z=freq['z'])
            stats = get_predictive_stats(samples=predictive_samples)
            stats.update(get_discrete_stats(samples=discrete_samples))
            output_by_type[svtype] = get_output(vids_np=vids_np, genotypes=genotypes, stats=stats, args=args)
    results_output = {}
    for svtype in output_by_type:
        results_output.update(output_by_type[svtype])
    return results_output


def get_output(vids_np: np.ndarray, genotypes: dict, stats: dict, args):
    n_variants = vids_np.size
    output_dict = {}
    for i in range(n_variants):
        vid = vids_np[i]
        if 'eps_pe' in stats:
            eps_pe = stats['eps_pe']['mean'][i]
        else:
            eps_pe = 0
        output_dict[vid] = {
            'gt': genotypes['gt'][i, :],
            'gt_p': genotypes['gt_p'][i, :],
            'gt_lod': genotypes['gt_lod'][i, :],
            'p_m_pe': stats['m_pe']['mean'][i],
            'p_m_sr1': stats['m_sr1']['mean'][i],
            'p_m_sr2': stats['m_sr2']['mean'][i],
            'p_m_rd': stats['m_rd']['mean'][i],
            'eps_pe': eps_pe * args.eps_pe,
            'eps_sr1': stats['eps_sr1']['mean'][i] * args.eps_sr1,
            'eps_sr2': stats['eps_sr2']['mean'][i] * args.eps_sr2
        }
    return output_dict


def get_global_stats(stats: dict):
    return {
        'pi_pe_mean': stats['pi_pe']['mean'],
        'pi_pe_std': stats['pi_pe']['std'],
        'pi_sr1_mean': stats['pi_sr1']['mean'],
        'pi_sr1_std': stats['pi_sr1']['std'],
        'pi_sr2_mean': stats['pi_sr2']['mean'],
        'pi_sr2_std': stats['pi_sr2']['std'],
        'pi_rd_mean': stats['pi_rd']['mean'],
        'pi_rd_std': stats['pi_rd']['std']
    }


def get_predictive_stats(samples: dict):
    return {key: {'mean': samples[key].astype(dtype='float').mean(axis=0).squeeze()} for key in samples}


def get_discrete_stats(samples: dict):
    return {key: {'mean': samples[key].astype(dtype='float').mean(axis=0).squeeze(),
                  'std': samples[key].astype(dtype='float').std(axis=0).squeeze()} for key in ['m_pe', 'm_sr1', 'm_sr2', 'm_rd']}


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
    z_bins = np.arange(model.k + 1) - 0.5
    for i in range(z.shape[1]):
        for j in range(z.shape[2]):
            z_freq[i, j, :], _ = np.histogram(z[:, i, j], bins=z_bins, density=True)
    m_pe_freq = discrete_samples['m_pe'].mean(axis=0)
    m_sr1_freq = discrete_samples['m_sr1'].mean(axis=0)
    m_sr2_freq = discrete_samples['m_sr2'].mean(axis=0)
    return {
        'z': z_freq,
        'm_pe': m_pe_freq,
        'm_sr1': m_sr1_freq,
        'm_sr2': m_sr2_freq
    }
