package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.List;

public class GenotypingEngine {
    /**
     * Function that fills vector with allele frequency priors. By default, infinite-sites, neutral variation prior is used,
     * where Pr(AC=i) = theta/i where theta is heterozygosity
     * @param N                                Number of chromosomes
     * @param priors                           (output) array to be filled with priors
     * @param heterozygosity                   default heterozygosity to use, if inputPriors is empty
     * @param inputPriors                      Input priors to use (in which case heterozygosity is ignored)
     */
    public static void computeAlleleFrequencyPriors(final int N, final double[] priors, final double heterozygosity, final List<Double> inputPriors) {


        double sum = 0.0;

        if (!inputPriors.isEmpty()) {
            // user-specified priors
            if (inputPriors.size() != N) {
                throw new UserException.BadArgumentValue("inputPrior", "Invalid length of inputPrior vector: vector length must be equal to # samples +1 ");
            }

            int idx = 1;
            for (final double prior: inputPriors) {
                if (prior < 0.0) {
                    throw new UserException.BadArgumentValue("Bad argument: negative values not allowed", "inputPrior");
                }
                priors[idx++] = Math.log10(prior);
                sum += prior;
            }
        }
        else {
            // for each i
            for (int i = 1; i <= N; i++) {
                final double value = heterozygosity / (double)i;
                priors[i] = Math.log10(value);
                sum += value;
            }
        }

        // protection against the case of heterozygosity too high or an excessive number of samples (which break population genetics assumptions)
        if (sum > 1.0) {
            throw new UserException.BadArgumentValue("heterozygosity","The heterozygosity value is set too high relative to the number of samples to be processed, or invalid values specified if input priors were provided - try reducing heterozygosity value or correct input priors.");
        }
        // null frequency for AF=0 is (1 - sum(all other frequencies))
        priors[0] = Math.log10(1.0 - sum);
    }

}
