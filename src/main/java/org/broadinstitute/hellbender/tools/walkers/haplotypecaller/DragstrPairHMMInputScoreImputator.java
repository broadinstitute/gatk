package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrReadSTRAnalyzer;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Pair-HMM score imputator based on the DRAGstr model parameters.
 * <p>
 *     The match, gap-openning and gap-extension at each position is chosen
 *     based on the sequence context (STR unit and repeat length) around that base on the read-sequence.
 * </p>
 */
public class DragstrPairHMMInputScoreImputator implements PairHMMInputScoreImputator {

    /**
     * As per the matlab scripts provided the GOP and GCP assigned to the last position of the read are 45 and 10
     * respectively.
     */
    private static final int GOP_AT_THE_END_OF_READ = 45;
    private static final int GCP_AT_THE_END_OF_READ = 10;

    /**
     * Cap for GOP values within the read.
     * <p>
     *   Therefore in practice GOP larger than this number are capped to its value.
     * </p>
     */
    private static final double MAX_GOP_IN_READ = 40.0;

    /**
     * Holds a reference to the Dragstr Model parameters.
     */
    private final DragstrParams params;

    /**
     * Constructs an inputator given the DragstrParameter collection.
     * @param params params dragstr model parameter collection.
     * @throws org.broadinstitute.hellbender.exceptions.UserException.CouldNotReadInputFile if there was a problem
     * accessing the parameter file.
     * @return never {@code null}.
     */
    public static DragstrPairHMMInputScoreImputator of(final DragstrParams params) {
        return new DragstrPairHMMInputScoreImputator(params);
    }

    private DragstrPairHMMInputScoreImputator(final DragstrParams params) {
        this.params = params;
    }

    @Override
    public PairHMMInputScoreImputation impute(final GATKRead read) {
        final byte[] bases = read.getBases();
        final int length = bases.length;
        final byte[] gop = new byte[length];
        final byte[] gcp = new byte[length];

        final DragstrReadSTRAnalyzer analyzer =  DragstrReadSTRAnalyzer.of(read, params.maximumPeriod());

        for (int i = 0; i < length - 1; i++) {
            final int period = analyzer.mostRepeatedPeriod(i);
            final int repeats = analyzer.numberOfMostRepeats(i);
            gop[i] = (byte) Math.round(Math.min(MAX_GOP_IN_READ, params.gop(period, repeats)));
            gcp[i] = (byte) Math.round(params.gcp(period, repeats));
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
