package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.pairhmm.DragstrParams;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class DragstrPairHMMInputScoreImputator implements PairHMMInputScoreImputator {

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
        final DragstrUtils.STRSequenceAnalyzer analyzer = DragstrUtils.repeatPeriodAndCounts(read.getLength(), params.maximumPeriod());
        analyzer.load(bases);
        final int length = bases.length;
        final byte[] gop = new byte[length];
        final byte[] gcp = new byte[length];
        for (int i = 0; i < length; i++) {
            final int period = analyzer.mostRepeatedPeriod(i);
            final int repeats = analyzer.numberOfMostRepeats(i);
            gop[i] = (byte) params.gop(period, repeats);
            gcp[i] = (byte) params.gcp(period, repeats);
        }
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
