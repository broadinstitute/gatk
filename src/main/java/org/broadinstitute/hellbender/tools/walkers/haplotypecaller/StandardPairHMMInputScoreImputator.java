package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Arrays;

/**
 * Standard or classic pair-hmm score imputator.
 * <p>
 *     This implements the default score calculations before the introduction of DRAGstr.
 * </p>
 */
public class StandardPairHMMInputScoreImputator implements PairHMMInputScoreImputator {

    private final byte constantGCP;

    private StandardPairHMMInputScoreImputator(final byte constantGCP) {
        this.constantGCP = constantGCP;

    }

    public static StandardPairHMMInputScoreImputator newInstance(final byte constantGCP) {
        return new StandardPairHMMInputScoreImputator(constantGCP);
    }

    @Override
    public PairHMMInputScoreImputation impute(final GATKRead read) {
        return new PairHMMInputScoreImputation() {
            @Override
            public byte[] delOpenPenalties() {
                return ReadUtils.getBaseDeletionQualities(read);
            }

            @Override
            public byte[] insOpenPenalties() {
                return ReadUtils.getBaseInsertionQualities(read);
            }

            @Override
            public byte[] gapContinuationPenalties() {
                final byte[] result = new byte[read.getLength()];
                Arrays.fill(result, constantGCP);
                return result;
            }
        };
    }
}