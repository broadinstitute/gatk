import numpy as np

import pyro
from pyro.infer import JitTraceEnum_ELBO, TraceEnum_ELBO
from pyro.infer.svi import SVI
from pyro.optim import PyroOptim

import torch
from torch.distributions import constraints

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
    pyro.clear_param_store()
    n_latent_dim = model.get_latent_dim(data.pe_t.shape[0])
    loc_init = torch.zeros(n_latent_dim, device=model.device)
    pyro.param("auto_loc", loc_init)
    scale_init = torch.ones(n_latent_dim, device=model.device)
    pyro.param("auto_scale", scale_init, constraint=constraints.positive)
    return svi.loss(model.model, model.guide, data_pe=data.pe_t, data_sr1=data.sr1_t, data_sr2=data.sr2_t,
                    depth_t=data.depth_t, rd_gt_prob_t=data.rd_gt_prob_t)


def run_training(model: SVGenotyperPyroModel,
                 data: SVGenotyperData,
                 scheduler: PyroOptim,
                 loss,
                 args) -> List[float]:
    logging.info("Running model training...")
    svi = SVI(model.model, model.guide, optim=scheduler, loss=loss)
    initialize_guide(svi=svi, model=model, data=data)
    train_elbo = []

    # Run training loop.  Use try to allow for keyboard interrupt.
    try:
        for epoch in range(args.max_iter):
            # Train, and keep track of training loss.
            total_epoch_loss = train_epoch(svi=svi, data=data, epoch=epoch, scheduler=scheduler)
            train_elbo.append(-total_epoch_loss)
            model.loss['train']['epoch'].append(epoch)
            model.loss['train']['elbo'].append(-total_epoch_loss)
            if epoch % args.iter_log_freq == 0:
                logging.info("[epoch %04d]  training loss: %.4f" % (epoch, total_epoch_loss))
        logging.info("Training procedure complete.")

    # Exception allows program to continue after ending training prematurely.
    except KeyboardInterrupt:
        logging.info("Training procedure stopped by keyboard interrupt.")

    return train_elbo


def run(args):
    logging.basicConfig(level=logging.INFO)
    pyro.enable_validation(True)
    pyro.distributions.enable_validation(True)
    pyro.clear_param_store()

    np.random.seed(args.random_seed)
    torch.random.manual_seed(args.random_seed)
    pyro.set_rng_seed(args.random_seed)

    scheduler = create_scheduler(lr_min=args.lr_min, lr_init=args.lr_init, lr_decay=args.lr_decay,
                                 beta1=args.adam_beta1, beta2=args.adam_beta2)

    if args.jit:
        elbo = JitTraceEnum_ELBO()
    else:
        elbo = TraceEnum_ELBO()

    for svtype in [SVTypes.DEL, SVTypes.DUP, SVTypes.INS, SVTypes.INV]:
        model = SVGenotyperPyroModel(svtype, args.num_states, args.eps_pe, args.eps_sr1, args.eps_sr2, args.device)
        data = load_data(vcf_path=args.vcf, mean_coverage_path=args.coverage_file, svtype=svtype, num_states=args.num_states,
                         depth_dilution_factor=args.depth_dilution_factor, device=args.device)
        logging.info("Training {:s} with {:d} variants and {:d} samples...".format(str(svtype), data.pe_t.shape[0], data.pe_t.shape[1]))
        train_elbo = run_training(model=model, data=data, scheduler=scheduler, loss=elbo, args=args)

