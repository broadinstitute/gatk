package org.broadinstitute.hellbender.tools.exome.pulldown;

import org.apache.commons.lang3.tuple.Pair;

import java.util.Collection;

/**
 * Balanced model prior for heterozygous pileups, i.e. allele fraction = 1/2.
 *
 * This prior is suitable for detecting heterozygous sites when the ref to alt allele ratio is expected to be 1/2,
 * which is the case for sequenced reads from pure germline samples.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class BalancedHeterozygousPileupPriorModel extends HeterozygousPileupPriorModel {

    /**
     * Calculates the log likelihood of heterozygosity assuming allele fraction = 1/2
     * @param coeffs list of (alpha, beta) tuples
     * @return any double value
     */
    @Override
    public double getHetLogLikelihood(final Collection<? extends Pair<Double, Double>> coeffs) {
        return getHetLogLikelihoodFixedAlleleFraction(0.5, coeffs);
    }
}
