import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.ops.indexing import Vindex
from pyro.infer import config_enumerate
from pyro.infer.autoguide import AutoDiagonalNormal

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
                 device: str = 'cpu'):
        self.k = k
        self.mu_eps_pe = mu_eps_pe
        self.mu_eps_sr1 = mu_eps_sr1
        self.mu_eps_sr2 = mu_eps_sr2
        self.svtype = svtype
        self.device = device
        self.loss = {'train': {'epoch': [], 'elbo': []},
                     'test': {'epoch': [], 'elbo': []}}

        if svtype == SVTypes.DEL or svtype == SVTypes.DUP:
            self.latent_sites = ['eps_pe', 'eps_sr1', 'eps_sr2', 'pi_pe', 'pi_sr1', 'pi_sr2', 'pi_rd']
        elif svtype == SVTypes.INS:
            self.latent_sites = ['eps_sr1', 'eps_sr2', 'pi_sr1', 'pi_sr2']
        elif svtype == SVTypes.INV:
            self.latent_sites = ['eps_pe', 'eps_sr1', 'eps_sr2', 'pi_sr1', 'pi_sr2', 'pi_pe']
        else:
            raise ValueError('SV type {:s} not supported for genotyping.'.format(str(svtype)))

        self.guide = AutoDiagonalNormal(poutine.block(self.model, expose=self.latent_sites))

    def get_latent_dim(self, n_variants: int):
        if self.svtype == SVTypes.DEL or self.svtype == SVTypes.DUP:
            return 4 + n_variants * 3
        elif self.svtype == SVTypes.INS:
            return 2 + n_variants * 2
        elif self.svtype == SVTypes.INV:
            return 3 + n_variants * 3
        else:
            raise ValueError('SV type {:s} not supported for genotyping.'.format(str(self.svtype)))

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

        if self.svtype == SVTypes.DEL or self.svtype == SVTypes.DUP:
            pi_rd = pyro.sample('pi_rd', dist.Beta(one_t, one_t))

        if self.svtype != SVTypes.INS:
            pi_pe = pyro.sample('pi_pe', dist.Beta(one_t, one_t))

        with pyro.plate('variant', n_variants, dim=-2, device=self.device):
            if self.svtype == SVTypes.DEL or self.svtype == SVTypes.DUP:
                m_rd = pyro.sample('m_rd', dist.Bernoulli(pi_rd))
            else:
                m_rd = zero_t.expand(n_variants)

            if self.svtype != SVTypes.INS:
                eps_pe = pyro.sample('eps_pe', dist.Exponential(one_t)).unsqueeze(-1) * self.mu_eps_pe
                eps_expanded_pe = eps_pe.expand(n_variants, n_samples, 1)
                m_pe = pyro.sample('m_pe', dist.Bernoulli(pi_pe))
            else:
                m_pe = zero_t.expand(n_variants)

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
                    # V x S
                    pyro.sample('obs_pe', dist.Poisson(rate=mu_obs_pe), obs=data_pe)
                else:
                    mu_obs_pe = zero_t.expand(n_variants, n_samples)

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
                # V x S
                pyro.sample('obs_sr1', dist.Poisson(rate=mu_obs_sr1), obs=data_sr1)
                pyro.sample('obs_sr2', dist.Poisson(rate=mu_obs_sr2), obs=data_sr2)
        return {
            'z': z,
            'm_pe': m_pe,
            'm_sr1': m_sr1,
            'm_sr2': m_sr2,
            'm_rd': m_rd
        }
