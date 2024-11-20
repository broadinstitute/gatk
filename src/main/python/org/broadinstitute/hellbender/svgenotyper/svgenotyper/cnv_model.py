import logging
import numpy as np
import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.ops.indexing import Vindex
from pyro.infer import config_enumerate, infer_discrete
from pyro.infer.predictive import Predictive
from pyro.infer.autoguide import AutoDiagonalNormal
import torch


class SVDepthData(object):
    def __init__(self,
                 samples: np.ndarray,
                 sample_per_base_depth: torch.Tensor,
                 sample_ploidy: torch.Tensor,
                 contigs: np.ndarray,
                 starts: torch.Tensor,
                 bin_size: torch.Tensor,
                 counts: torch.Tensor):
        self.samples = samples
        self.sample_per_base_depth = sample_per_base_depth
        self.sample_ploidy = sample_ploidy
        self.contigs = contigs
        self.starts = starts
        self.bin_size = bin_size
        self.counts = counts


class SVDepthPyroModel(object):
    def __init__(self,
                 k: int,
                 tensor_dtype: torch.dtype,
                 read_length: float,
                 var_phi_bin: float = 0.2,
                 var_phi_sample: float = 0.01,
                 alpha_non_ref: float = 1.,
                 alpha_ref: float = 18.,
                 mu_eps: float = 0.01,
                 device: str = 'cpu',
                 loss: dict = None):
        if (k != 5):  # TODO
            raise ValueError('Only k=5 is supported')
        self.k = k
        self.tensor_dtype = tensor_dtype
        self.var_phi_bin = var_phi_bin
        self.var_phi_sample = var_phi_sample
        self.mu_eps = mu_eps
        self.alpha_non_ref = alpha_non_ref
        self.alpha_ref = alpha_ref
        self.read_length = read_length
        self.device = device
        if loss is None:
            self.loss = {'epoch': [], 'elbo': []}
        else:
            self.loss = loss

        self.latent_sites = ['phi_bin', 'phi_sample', 'eps', 'p_hw']
        self.guide = AutoDiagonalNormal(poutine.block(self.model, expose=self.latent_sites))

    @config_enumerate(default="parallel")
    def model(self,
              counts: torch.Tensor,
              bin_size: torch.Tensor,
              sample_ploidy: torch.Tensor,
              sample_depth: torch.Tensor):
        n_intervals = counts.shape[0]
        n_samples = counts.shape[1]

        if sample_ploidy[sample_ploidy > 2].sum() > 0:
            raise ValueError("Sample ploidy >2 not supported")

        zero_t = torch.zeros(1, device=self.device, dtype=self.tensor_dtype)
        one_t = torch.ones(1, device=self.device, dtype=self.tensor_dtype)
        k_range_t = torch.arange(0, self.k).to(device=self.device, dtype=self.tensor_dtype)

        plate_1 = pyro.plate('interval', n_intervals, dim=-2, device=self.device)
        plate_2 = pyro.plate('sample', n_samples, dim=-1, device=self.device)

        with plate_2:
            phi_sample = pyro.sample('phi_sample', dist.LogNormal(zero_t, self.var_phi_sample))

        with plate_1:

            phi_bin = pyro.sample('phi_bin', dist.LogNormal(zero_t, self.var_phi_bin))
            eps_rd = pyro.sample('eps', dist.Exponential(one_t)) * self.mu_eps

            alpha_hw = self.alpha_non_ref * one_t.expand(3)
            alpha_hw[1] = self.alpha_ref
            p_hw = pyro.sample('p_hw', dist.Dirichlet(alpha_hw))

        with plate_1, plate_2:
            locs = phi_bin.unsqueeze(-1) * k_range_t.unsqueeze(0).unsqueeze(0) + eps_rd.unsqueeze(-1)

            p_hw_loss = p_hw[..., 0]
            p_hw_neutral = p_hw[..., 1]
            p_hw_gain = p_hw[..., 2]
            z0 = p_hw_loss * p_hw_loss
            z1 = 2 * p_hw_loss * p_hw_neutral
            z2 = p_hw_neutral * p_hw_neutral + 2 * p_hw_loss * p_hw_gain
            z3 = 2 * p_hw_gain * p_hw_neutral
            z4 = p_hw_gain * p_hw_gain

            z0_hap = p_hw_loss
            z1_hap = p_hw_neutral
            z2_hap = p_hw_gain
            z3_hap = zero_t.unsqueeze(-1).expand(n_intervals, 1)
            z4_hap = zero_t.unsqueeze(-1).expand(n_intervals, 1)

            z_dist_dip = torch.stack([z0, z1, z2, z3, z4], dim=-1).expand(n_intervals, n_samples, self.k)
            z_dist_hap = torch.stack([z0_hap, z1_hap, z2_hap, z3_hap, z4_hap], dim=-1).expand(n_intervals, n_samples, self.k)
            # Apply haploid model for ploidy 0 or 1, otherwise diploid
            z_dist = torch.where(sample_ploidy.unsqueeze(-1) < 2, z_dist_hap, z_dist_dip)
            z = pyro.sample('z', dist.Categorical(z_dist))

            # N x S
            # TODO add in * phi_sample
            mu_counts = bin_size.unsqueeze(-1) * sample_depth * Vindex(locs)[..., z] / self.read_length

            pyro.sample('obs', dist.Poisson(rate=mu_counts), obs=counts)

    def run_predictive(self, data: SVDepthData, n_samples: int = 100, n_iter: int = 10):
        logging.info("Running predictive distribution inference...")
        predictive = Predictive(self.model, guide=self.guide, num_samples=n_samples, return_sites=self.latent_sites)

        posterior_means = None
        for i in range(n_iter):
            logging.info("Iteration {:d} of {:d}...".format(i + 1, n_iter))
            sample = predictive(counts=data.counts, bin_size=data.bin_size, sample_ploidy=data.sample_ploidy,
                                sample_depth=data.sample_per_base_depth)
            sample_means = {key: {"mean": sample[key].detach().cpu().numpy().astype(dtype='float').mean(axis=0).squeeze()} for key in sample}
            if posterior_means is None:
                posterior_means = sample_means
            else:
                for key in posterior_means:
                    posterior_means[key]["mean"] += sample_means[key]["mean"]

        for key in posterior_means:
            posterior_means[key]["mean"] = posterior_means[key]["mean"] / n_iter
        logging.info("Inference complete.")
        return posterior_means

    def run_discrete(self, data: SVDepthData, n_states: int, log_freq: int = 100, n_samples: int = 1000):
        logging.info("Running discrete inference...")
        posterior_means = None
        guide_trace = poutine.trace(self.guide).get_trace(counts=data.counts, bin_size=data.bin_size, sample_ploidy=data.sample_ploidy,
                                                          sample_depth=data.sample_per_base_depth)
        trained_model = poutine.replay(self.model, trace=guide_trace)
        with torch.no_grad():
            for i in range(n_samples):
                inferred_model = infer_discrete(trained_model, temperature=1, first_available_dim=-3)
                trace = poutine.trace(inferred_model).get_trace(counts=data.counts, bin_size=data.bin_size,
                                                                sample_ploidy=data.sample_ploidy,
                                                                sample_depth=data.sample_per_base_depth)
                states_sample = trace.nodes["z"]["value"].detach().cpu().numpy().astype(dtype='float').squeeze()
                states_sample_one_hot = np.zeros([states_sample.shape[0], states_sample.shape[1], n_states], dtype='float')
                for j in range(n_states):
                    states_sample_one_hot[states_sample == j, j] = 1.
                if posterior_means is None:
                    posterior_means = {"z": {"mean": states_sample_one_hot}}
                else:
                    posterior_means["z"]["mean"] += states_sample_one_hot
                if (i + 1) % log_freq == 0:
                    logging.info("[sample {:d}] discrete latent".format(i + 1))
        for key in posterior_means:
            posterior_means[key]["mean"] = posterior_means[key]["mean"] / n_samples
        logging.info("Inference complete.")
        return posterior_means
