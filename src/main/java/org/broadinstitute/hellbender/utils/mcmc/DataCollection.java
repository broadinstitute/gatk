package org.broadinstitute.hellbender.utils.mcmc;

/**
 * Interface for tagging any class that represents a collection of datasets required to update posterior samples
 * for Markov-Chain Monte Carlo sampling using samplers implementing the {@link ParameterSampler} interface.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public interface DataCollection {}
