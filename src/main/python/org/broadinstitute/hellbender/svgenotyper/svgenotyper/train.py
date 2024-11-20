import json
import logging
import numpy as np
import os
import pyro
from pyro.infer import JitTraceEnum_ELBO, TraceEnum_ELBO
from pyro.infer.svi import SVI
from pyro.optim import PyroOptim
import torch

from .constants import SVTypes
from .model import SVGenotyperData, SVGenotyperPyroModel
from . import io


def train_epoch(svi: SVI, data: SVGenotyperData, epoch: int, scheduler: pyro.optim.LambdaLR) -> float:
    loss = svi.step(data_pe=data.pe_t, data_sr1=data.sr1_t, data_sr2=data.sr2_t, depth_t=data.depth_t,
                    svlen_t=data.svlen_t, rd_gt_prob_t=data.rd_gt_prob_t)
    scheduler.step(epoch)
    return loss


def create_scheduler(lr_min: float, lr_init: float, lr_decay: float, beta1: float, beta2: float):
    return pyro.optim.LambdaLR({
        'optimizer': torch.optim.Adam,
        'optim_args': {'lr': 1., 'betas': (beta1, beta2)},
        'lr_lambda': lambda k: lr_min + (lr_init - lr_min) * np.exp(-k / lr_decay)
    })


def run_training(model: SVGenotyperPyroModel,
                 data: SVGenotyperData,
                 max_iter: int,
                 lr_min: float,
                 lr_init: float,
                 lr_decay: float,
                 adam_beta1: float,
                 adam_beta2: float,
                 iter_log_freq: int,
                 jit: bool):
    logging.info("Initializing model...")
    pyro.clear_param_store()
    scheduler = create_scheduler(lr_min=lr_min, lr_init=lr_init, lr_decay=lr_decay,
                                 beta1=adam_beta1, beta2=adam_beta2)
    if jit:
        loss = JitTraceEnum_ELBO()
    else:
        loss = TraceEnum_ELBO()
    svi = SVI(model.model, model.guide, optim=scheduler, loss=loss)
    logging.info("Running model training...")
    # Run training loop.  Use try to allow for keyboard interrupt.
    try:
        for epoch in range(max_iter):
            # Train, and keep track of training loss.
            total_epoch_loss = train_epoch(svi=svi, data=data, epoch=epoch, scheduler=scheduler)
            model.loss['epoch'].append(epoch)
            model.loss['elbo'].append(-total_epoch_loss)
            if (epoch + 1) % iter_log_freq == 0:
                logging.info("[epoch %04d]  training loss: %.4f" % (epoch + 1, total_epoch_loss))
        logging.info("Training procedure complete after {:d} epochs, loss: {:.4f}".format(max_iter, total_epoch_loss))

    # Exception allows program to continue after ending training prematurely.
    except KeyboardInterrupt:
        logging.info("Training procedure stopped by keyboard interrupt.")


def run(args: dict, default_dtype: torch.dtype = torch.float32):
    base_path = os.path.join(args['output_dir'], args['output_name'])
    log_path = base_path + ".log.txt"
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

    model = SVGenotyperPyroModel(svtype=svtype, k=args['num_states'], tensor_dtype=default_dtype, mu_eps_pe=args['eps_pe'], mu_eps_sr1=args['eps_sr1'],
                                 mu_eps_sr2=args['eps_sr2'], var_phi_pe=args['phi_pe'], var_phi_sr1=args['phi_sr1'],
                                 var_phi_sr2=args['phi_sr2'], read_length=args['read_length'], device=args['device'])
    data = io.load_data(variants_file_path=args['variants_file'], mean_coverage_path=args['coverage_file'], samples_path=args['samples_file'],
                        svtype=svtype, num_states=model.k, depth_dilution_factor=args['depth_dilution_factor'],
                        tensor_dtype=default_dtype, device=args['device'])
    if data is None:
        logging.warning("No variants found.")
    else:
        logging.info("Training {:s} with {:d} variants and {:d} samples...".format(str(svtype.name), data.pe_t.shape[0], data.pe_t.shape[1]))
        run_training(model=model, data=data, max_iter=args['max_iter'], lr_min=args['lr_min'], lr_init=args['lr_init'],
                     lr_decay=args['lr_decay'], adam_beta1=args['adam_beta1'], adam_beta2=args['adam_beta2'],
                     iter_log_freq=args['iter_log_freq'], jit=args['jit'])
        save_data(base_path=base_path, model=model, data=data, vids=data.vids, sample_ids=data.samples)


def save_data(base_path: str, model: SVGenotyperPyroModel, data: SVGenotyperData, vids: list, sample_ids: list):
    io.save_tensors(data=data, base_path=base_path)
    io.save_list(data=vids, path=base_path + ".vids.list")
    io.save_list(data=sample_ids, path=base_path + ".sample_ids.list")
    save_model(model=model, base_path=base_path)


def save_model(model: SVGenotyperPyroModel, base_path: str):
    with open(base_path + '.vars.json', 'w') as f:
        json.dump({
            'k': model.k,
            'mu_eps_pe': model.mu_eps_pe,
            'mu_eps_sr1': model.mu_eps_sr1,
            'mu_eps_sr2': model.mu_eps_sr2,
            'var_phi_pe': model.var_phi_pe,
            'var_phi_sr1': model.var_phi_sr1,
            'var_phi_sr2': model.var_phi_sr2,
            'svtype': model.svtype.name,
            'read_length': model.read_length,
            'loss': model.loss
        }, f)
    pyro.get_param_store().save(base_path + '.param_store.pyro')
