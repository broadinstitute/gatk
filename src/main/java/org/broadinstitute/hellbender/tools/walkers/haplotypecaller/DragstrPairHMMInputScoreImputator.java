package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.pairhmm.DragstrParams;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrReadSTRAnalyzer;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Pair-HMM score imputator based on the DRAGstr model parameters.
 * <p>
 *     The match, gap-openning and gap-extension at each position is chosen
 *     based on the sequence context (STR unit and repeat length) around that base on the read-sequence.
 * </p>
 */
public class DragstrPairHMMInputScoreImputator implements PairHMMInputScoreImputator {

    private static final int GOP_AT_THE_END_OF_READ = 45;
    private static final int GCP_AT_THE_END_OF_READ = 10;

    private final DragstrParams params;

    public static DragstrPairHMMInputScoreImputator newInstance(final String path) {
        return new DragstrPairHMMInputScoreImputator(new DragstrParams(path));
    }

    public static DragstrPairHMMInputScoreImputator newInstance(final DragstrParams params) {
        return new DragstrPairHMMInputScoreImputator(params);
    }

    private DragstrPairHMMInputScoreImputator(final DragstrParams params) {
        this.params = params;
    }

    @Override
    public PairHMMInputScoreImputation impute(final GATKRead read) {
        final byte[] bases = read.getBases();
        final DragstrReadSTRAnalyzer analyzer = DragstrUtils.repeatPeriodAndCounts(read.getLength(),
                params.maximumPeriod());
        analyzer.load(bases);
        final int length = bases.length;
        final byte[] gop = new byte[length];
        final byte[] gcp = new byte[length];
        for (int i = 0; i < length - 1; i++) {
            final int period = analyzer.mostRepeatedPeriod(i + 1);
            final int repeats = analyzer.numberOfMostRepeats(i + 1);
            gop[i] = (byte) params.gop(period, repeats);
            gcp[i] = (byte) params.gcp(period, repeats);
        }
        gop[length - 1] = GOP_AT_THE_END_OF_READ;
        gcp[length - 1] = GCP_AT_THE_END_OF_READ;
        return new PairHMMInputScoreImputation() {
            @Override
            public byte[] delOpenPenalties() {
                return gop;
            }

            @Override
            public byte[] insOpenPenalties() {
                return gop;
            }

            @Override
            public byte[] gapContinuationPenalties() {
                return gcp;
            }
        };
    }
}
