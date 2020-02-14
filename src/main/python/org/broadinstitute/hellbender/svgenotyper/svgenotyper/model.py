import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.ops.indexing import Vindex
from pyro.infer import config_enumerate, Predictive, infer_discrete
from pyro.infer.autoguide import AutoDiagonalNormal

import logging
import torch

from .constants import SVTypes


class SVGenotyperData(object):
    def __init__(self,
                 pe_t: torch.Tensor,
                 sr1_t: torch.Tensor,
                 sr2_t: torch.Tensor,
                 depth_t: torch.Tensor,
                 rd_gt_prob_t: torch.Tensor):
        self.pe_t = pe_t
        self.sr1_t = sr1_t
        self.sr2_t = sr2_t
        self.depth_t = depth_t
        self.rd_gt_prob_t = rd_gt_prob_t


class SVGenotyperPyroModel(object):
    def __init__(self,
                 svtype: SVTypes,
                 k: int = 5,
                 mu_eps_pe: float = 0.1,
                 mu_eps_sr1: float = 0.1,
                 mu_eps_sr2: float = 0.1,
                 mu_lambda_pe: float = 0.1,
                 mu_lambda_sr1: float = 0.1,
                 mu_lambda_sr2: float = 0.1,
                 device: str = 'cpu'):
        self.k = k
        self.mu_eps_pe = mu_eps_pe
        self.mu_eps_sr1 = mu_eps_sr1
        self.mu_eps_sr2 = mu_eps_sr2
        self.mu_lambda_pe = mu_lambda_pe
        self.mu_lambda_sr1 = mu_lambda_sr1
        self.mu_lambda_sr2 = mu_lambda_sr2
        self.svtype = svtype
        self.device = device
        self.loss = {'train': {'epoch': [], 'elbo': []},
                     'test': {'epoch': [], 'elbo': []}}

        if svtype == SVTypes.DEL or svtype == SVTypes.DUP:
            self.latent_sites = ['eps_pe', 'eps_sr1', 'eps_sr2', 'pi_pe', 'pi_sr1', 'pi_sr2', 'pi_rd', 'var_pe',
                                 'var_sr1', 'var_sr2']
        elif svtype == SVTypes.INS:
            self.latent_sites = ['eps_sr1', 'eps_sr2', 'pi_sr1', 'pi_sr2', 'var_sr1', 'var_sr2']
        elif svtype == SVTypes.INV:
            self.latent_sites = ['eps_pe', 'eps_sr1', 'eps_sr2', 'pi_sr1', 'pi_sr2', 'pi_pe' 'var_pe', 'var_sr1',
                                 'var_sr2']
        else:
            raise ValueError('SV type {:s} not supported for genotyping.'.format(str(svtype.name)))

        self.guide = AutoDiagonalNormal(poutine.block(self.model, expose=self.latent_sites))

    def get_latent_dim(self, n_variants: int):
        if self.svtype == SVTypes.DEL or self.svtype == SVTypes.DUP:
            return 4 + n_variants * 3
        elif self.svtype == SVTypes.INS:
            return 2 + n_variants * 2
        elif self.svtype == SVTypes.INV:
            return 3 + n_variants * 3
        else:
            raise ValueError('SV type {:s} not supported for genotyping.'.format(str(self.svtype.name)))

    @config_enumerate(default="parallel")
    def model(self,
              data_pe: torch.Tensor,
              data_sr1: torch.Tensor,
              data_sr2: torch.Tensor,
              depth_t: torch.Tensor,
              rd_gt_prob_t: torch.Tensor):

        n_variants = data_pe.shape[0]
        n_samples = data_pe.shape[1]
        zero_t = torch.zeros(1, device=self.device)
        one_t = torch.ones(1, device=self.device)
        k_range_t = torch.arange(1, self.k).to(device=self.device)
        ones_vsk1_t = torch.ones(n_variants, n_samples, self.k - 1, device=self.device)
        w_uniform_t = torch.ones(self.k, device=self.device) / self.k

        pi_sr1 = pyro.sample('pi_sr1', dist.Beta(one_t, one_t))
        pi_sr2 = pyro.sample('pi_sr2', dist.Beta(one_t, one_t))

        var_sr1 = pyro.sample('var_sr1', dist.Exponential(one_t)) * self.mu_lambda_sr1
        var_sr2 = pyro.sample('var_sr2', dist.Exponential(one_t)) * self.mu_lambda_sr2

        if self.svtype == SVTypes.DEL or self.svtype == SVTypes.DUP:
            pi_rd = pyro.sample('pi_rd', dist.Beta(one_t, one_t))

        if self.svtype != SVTypes.INS:
            pi_pe = pyro.sample('pi_pe', dist.Beta(one_t, one_t))
            var_pe = pyro.sample('var_pe', dist.Exponential(one_t)) * self.mu_lambda_pe

        with pyro.plate('variant', n_variants, dim=-2, device=self.device):
            if self.svtype == SVTypes.DEL or self.svtype == SVTypes.DUP:
                m_rd = pyro.sample('m_rd', dist.Bernoulli(pi_rd))
            else:
                m_rd = zero_t.expand(n_variants).expand(n_variants).unsqueeze(-1)

            if self.svtype != SVTypes.INS:
                eps_pe = pyro.sample('eps_pe', dist.Exponential(one_t)).unsqueeze(-1) * self.mu_eps_pe
                eps_expanded_pe = eps_pe.expand(n_variants, n_samples, 1)
                m_pe = pyro.sample('m_pe', dist.Bernoulli(pi_pe))
            else:
                m_pe = zero_t.expand(n_variants).unsqueeze(-1)

            eps_sr1 = pyro.sample('eps_sr1', dist.Exponential(one_t)).unsqueeze(-1) * self.mu_eps_sr1
            eps_sr2 = pyro.sample('eps_sr2', dist.Exponential(one_t)).unsqueeze(-1) * self.mu_eps_sr2
            eps_expanded_sr1 = eps_sr1.expand(n_variants, n_samples, 1)
            eps_expanded_sr2 = eps_sr2.expand(n_variants, n_samples, 1)
            m_sr1 = pyro.sample('m_sr1', dist.Bernoulli(pi_sr1))
            m_sr2 = pyro.sample('m_sr2', dist.Bernoulli(pi_sr2))

            with pyro.plate('sample', n_samples, dim=-1, device=self.device):

                z_weights = m_rd.unsqueeze(-1) * rd_gt_prob_t + (one_t - m_rd.unsqueeze(-1)) * w_uniform_t
                z = pyro.sample('z', dist.Categorical(z_weights))

                if self.svtype != SVTypes.INS:
                    # V x S x (K-1)
                    nonzero_locs_pe = k_range_t * ones_vsk1_t
                    # V x S x K
                    m1_locs_pe = torch.cat([eps_expanded_pe, nonzero_locs_pe], dim=-1)
                    # V x S x K
                    m0_locs_pe = eps_expanded_pe.expand(n_variants, n_samples, self.k)
                    # V x S x K
                    gated_locs_pe = (one_t - m_pe.unsqueeze(-2)) * m0_locs_pe + m_pe.unsqueeze(-2) * m1_locs_pe
                    # V x S
                    mu_obs_pe = depth_t * Vindex(gated_locs_pe)[..., z]
                    var_pe = mu_obs_pe * (1. + var_pe)
                    r_pe = mu_obs_pe * mu_obs_pe / (var_pe - mu_obs_pe)
                    p_pe = (var_pe - mu_obs_pe) / var_pe
                    pyro.sample('obs_pe', dist.NegativeBinomial(total_count=r_pe, probs=p_pe), obs=data_pe)
                else:
                    r_pe = one_t.unsqueeze(-1).expand(n_variants, n_samples)
                    p_pe = one_t.unsqueeze(-1).expand(n_variants, n_samples)

                # V x S x (K-1)
                nonzero_locs_sr = k_range_t * ones_vsk1_t
                # V x S x K
                m1_locs_sr1 = torch.cat([eps_expanded_sr1, nonzero_locs_sr], dim=-1)
                m1_locs_sr2 = torch.cat([eps_expanded_sr2, nonzero_locs_sr], dim=-1)
                # V x S x K
                m0_locs_sr1 = eps_expanded_sr1.expand(n_variants, n_samples, self.k)
                m0_locs_sr2 = eps_expanded_sr2.expand(n_variants, n_samples, self.k)
                # V x S x K
                gated_locs_sr1 = (1. - m_sr1.unsqueeze(-2)) * m0_locs_sr1 + m_sr1.unsqueeze(-2) * m1_locs_sr1
                gated_locs_sr2 = (1. - m_sr2.unsqueeze(-2)) * m0_locs_sr2 + m_sr2.unsqueeze(-2) * m1_locs_sr2
                # V x S
                mu_obs_sr1 = depth_t * Vindex(gated_locs_sr1)[..., z]
                mu_obs_sr2 = depth_t * Vindex(gated_locs_sr2)[..., z]
                var_sr1 = mu_obs_sr1 * (1. + var_sr1)
                var_sr2 = mu_obs_sr2 * (1. + var_sr2)
                r_sr1 = mu_obs_sr1 * mu_obs_sr1 / (var_sr1 - mu_obs_sr1)
                r_sr2 = mu_obs_sr2 * mu_obs_sr2 / (var_sr2 - mu_obs_sr2)
                p_sr1 = (var_sr1 - mu_obs_sr1) / var_sr1
                p_sr2 = (var_sr2 - mu_obs_sr2) / var_sr2
                pyro.sample('sr1_obs', dist.NegativeBinomial(total_count=r_sr1, probs=p_sr1), obs=data_sr1)
                pyro.sample('sr2_obs', dist.NegativeBinomial(total_count=r_sr2, probs=p_sr2), obs=data_sr2)
        return {
            'z': z,
            'r_pe': r_pe,
            'r_sr1': r_sr1,
            'r_sr2': r_sr2,
            'p_pe': p_pe,
            'p_sr1': p_sr1,
            'p_sr2': p_sr2,
            'm_pe': m_pe,
            'm_sr1': m_sr1,
            'm_sr2': m_sr2,
            'm_rd': m_rd
        }

    def infer_predictive(self, data: SVGenotyperData, log_freq: int = 100, n_samples: int = 1000):
        logging.info("Running predictive distribution inference...")
        samples = []
        for i in range(n_samples):
            predictive = Predictive(self.model, guide=self.guide, num_samples=1, return_sites=self.latent_sites)
            sample = predictive(data_pe=data.pe_t, data_sr1=data.sr1_t, data_sr2=data.sr2_t, depth_t=data.depth_t, rd_gt_prob_t=data.rd_gt_prob_t)
            sample = {key: sample[key].detach().cpu() for key in sample}
            samples.append(sample)
            if (i + 1) % log_freq == 0:
                logging.info("[sample {:d}] predictive".format(i + 1))
        result = {key: torch.stack([samples[i][key] for i in range(n_samples)], dim=0).numpy() for key in self.latent_sites}
        logging.info("Inference complete.")
        return result

    def infer_discrete(self, data: SVGenotyperData, svtype: SVTypes, log_freq: int = 100, n_samples: int = 1000):
        logging.info("Running discrete inference...")
        posterior_samples = []
        for i in range(n_samples):
            guide_trace = poutine.trace(self.guide).get_trace(data_pe=data.pe_t, data_sr1=data.sr1_t,
                                                              data_sr2=data.sr2_t, depth_t=data.depth_t,
                                                              rd_gt_prob_t=data.rd_gt_prob_t)
            trained_model = poutine.replay(self.model, trace=guide_trace)
            inferred_model = infer_discrete(trained_model, temperature=1, first_available_dim=-3)
            trace = poutine.trace(inferred_model).get_trace(data_pe=data.pe_t, data_sr1=data.sr1_t,
                                                            data_sr2=data.sr2_t, depth_t=data.depth_t,
                                                            rd_gt_prob_t=data.rd_gt_prob_t)
            posterior_samples.append([trace.nodes["z"]["value"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["r_pe"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["r_sr1"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["r_sr2"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["p_pe"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["p_sr1"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["p_sr2"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["m_pe"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["m_sr1"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["m_sr2"].detach().cpu(),
                                      trace.nodes["_RETURN"]["value"]["m_rd"].detach().cpu()])
            if (i + 1) % log_freq == 0:
                logging.info("[sample {:d}] discrete latent".format(i + 1))
        posterior_samples = [torch.stack([posterior_samples[j][i] for j in range(n_samples)], dim=0) for i in range(11)]

        z = posterior_samples[0]
        r_pe = posterior_samples[1]
        r_sr1 = posterior_samples[2]
        r_sr2 = posterior_samples[3]
        p_pe = posterior_samples[4]
        p_sr1 = posterior_samples[5]
        p_sr2 = posterior_samples[6]
        m_pe = posterior_samples[7]
        m_sr1 = posterior_samples[8]
        m_sr2 = posterior_samples[9]
        m_rd = posterior_samples[10]

        samples_pe = []
        samples_sr1 = []
        samples_sr2 = []
        for i in range(n_samples):
            if svtype == SVTypes.INS:
                samples_pe.append(torch.zeros(1, device='cpu').unsqueeze(-1).expand(r_pe.shape[0], r_pe.shape[1]))
            else:
                samples_pe.append(dist.NegativeBinomial(total_count=r_pe[..., i], probs=p_pe[..., i]).sample())
            samples_sr1.append(dist.NegativeBinomial(total_count=r_sr1[..., i], probs=p_sr1[..., i]).sample())
            samples_sr2.append(dist.NegativeBinomial(total_count=r_sr2[..., i], probs=p_sr2[..., i]).sample())
            if (i + 1) % log_freq == 0:
                logging.info("[sample {:d}] discrete observed".format(i + 1))

        pe = torch.stack(samples_pe, dim=0)
        sr1 = torch.stack(samples_sr1, dim=0)
        sr2 = torch.stack(samples_sr2, dim=0)
        logging.info("Inference complete.")
        return {
            "z": z.numpy(),
            "pe": pe.numpy(),
            "sr1": sr1.numpy(),
            "sr2": sr2.numpy(),
            "m_pe": m_pe.numpy(),
            "m_sr1": m_sr1.numpy(),
            "m_sr2": m_sr2.numpy(),
            "m_rd": m_rd.numpy()
        }

