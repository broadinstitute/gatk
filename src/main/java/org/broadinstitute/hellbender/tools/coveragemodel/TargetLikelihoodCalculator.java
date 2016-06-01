package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.tools.exome.Target;

/**
 * Calculate log-likelihoods of observed coverage as a function of copy ratio.  This is needed for
 * likelihood-based segmentation and calling of CNVs, for example our HMM for germline CNVs.
 *
 * This interface is general in that "coverage" can represent any double associated with a target.
 * In particular, coverage could be "raw" integer read counts (cast to double), "proportional coverage"
 * normalized by average depth, log-coverage, or denoised "tangent-normalized" coverage obtained by
 * subtracting principal components.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@FunctionalInterface
public interface TargetLikelihoodCalculator {
    double logLikelihood(final Target target, final double copyRatio, final double coverage);
}
