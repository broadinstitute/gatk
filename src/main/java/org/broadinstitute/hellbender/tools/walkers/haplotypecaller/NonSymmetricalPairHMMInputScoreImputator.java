package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;

/**
 * A version of the classic {@link StandardPairHMMInputScoreImputator} that allows for decoupled insertion and deletion penalties for the model.
 * <p>
 *     This implements the default score calculations before the introduction of DRAGstr.
 *     This allows for different treatment of insertion and deletion penalties as flat scores
 * </p>
 */
public class NonSymmetricalPairHMMInputScoreImputator implements PairHMMInputScoreImputator {

    private final byte constantGCP;
    private final byte flatInsertionQual;
    private final byte flatDeletionQual;

    private NonSymmetricalPairHMMInputScoreImputator(final byte constantGCP, final byte flatInsertionQual, final byte flatDeletionQual) {
        this.constantGCP = constantGCP;
        this.flatDeletionQual = flatDeletionQual;
        this.flatInsertionQual = flatInsertionQual;

    }

    public static NonSymmetricalPairHMMInputScoreImputator newInstance(final byte constantGCP, final byte flatInsertionQual, final byte flatDeletionQual) {
        return new NonSymmetricalPairHMMInputScoreImputator(constantGCP, flatInsertionQual, flatDeletionQual);
    }

    @Override
    public PairHMMInputScoreImputation impute(final GATKRead read) {
        return new PairHMMInputScoreImputation() {
            @Override
            public byte[] delOpenPenalties() {
                byte[] quals = new byte[read.getBaseQualityCount()];
                Arrays.fill(quals, flatDeletionQual);  // Some day in the future when base insertion and base deletion quals exist the samtools API will
                // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
                return quals;
            }

            @Override
            public byte[] insOpenPenalties() {
                byte[] quals = new byte[read.getBaseQualityCount()];
                Arrays.fill(quals, flatInsertionQual);  // Some day in the future when base insertion and base deletion quals exist the samtools API will
                // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
                return quals;
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
