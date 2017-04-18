package org.broadinstitute.hellbender.tools.coveragemodel.interfaces;

import org.broadinstitute.hellbender.tools.exome.Target;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;

/**
 * Calculate log-likelihoods of observed coverage as a function of copy ratio.  This is needed for
 * likelihood-based segmentation and calling of CNVs, for example our HMM for germline CNVs.
 *
 * In this generic interface, {@code likelihoodData} can represent any type of data that is required
 * for calculating the likelihood function at a target. For example, the likelihood could represent
 * the emission probability to different copy number states.
 *
 * In particular, it could represent "raw" integer read counts (cast to double), "proportional coverage"
 * normalized by average depth, log-coverage, or denoised "tangent-normalized" coverage obtained by
 * subtracting noise principal components.
 *
 * @param <DATA> type of data required for calculating the likelihood
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@FunctionalInterface
public interface TargetLikelihoodCalculator<DATA> {
    double logLikelihood(@Nonnull final DATA datum, final double copyRatio, @Nullable final Target target);
}
