package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Arrays;

public class StandardPairHMMInputScoreImputator implements PairHMMInputScoreImputator {

    private final byte constantGCP;

    private final byte defaultGOP;

    private StandardPairHMMInputScoreImputator(final byte defaultGOP, final byte constantGCP) {
        this.constantGCP = constantGCP;
        this.defaultGOP = defaultGOP;
    }

    public static StandardPairHMMInputScoreImputator newInstance(final byte defaultGOP, final byte constantGCP) {
        return new StandardPairHMMInputScoreImputator(defaultGOP, constantGCP);
    }

    @Override
    public PairHMMInputScoreImputation impute(final GATKRead read) {
        return new PairHMMInputScoreImputation() {
            @Override
            public byte[] delOpenPenalties() {
                final byte[] existing = ReadUtils.getExistingBaseDeletionQualities(read);
                return existing == null ? defaultGOPs(read) : existing;
            }

            private byte[] defaultGOPs(final GATKRead read) {
                final byte[] existing = ReadUtils.getExistingBaseInsertionQualities(read);
                return existing == null ? defaultGOPs(read) : existing;
            }

            @Override
            public byte[] insOpenPenalties() {
                final byte[] existing = ReadUtils.getExistingBaseInsertionQualities(read);
                return existing == null ? defaultGOPs(read) : existing;
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
