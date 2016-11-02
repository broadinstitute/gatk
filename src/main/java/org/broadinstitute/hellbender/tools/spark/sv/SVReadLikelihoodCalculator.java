package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * An common read likelihood calculator interface for all types of SV's.
 *
 * All calculators are expected to be properly initialized, configured before run, and closed after running.
 *
 * Reads maybe preprocessed before sent for actual calculation; the default behavior is not to.
 */
interface SVReadLikelihoodCalculator {

    default void initialize() {}

    default void configure() {}

    default void close() {}

    /**
     * Scoring scheme for representing the probability of a read arising from a particular allele.
     */
    final class ScoringScheme implements Serializable{
        private static final long serialVersionUID = 1L;

        // TODO: these are the default setting for "-x intractg" as of bwa-mem 0.7.15-r1140
        public static final int MISMATCH_PEN_DEFAULT = 9;
        public static final int GAP_OPEN_PEN_DEFAULT = 16;
        public static final int GAP_EXTEN_PEN_DEFAULT = 1;
        public static final int CLIPPING_PEN_DEFAULT = 5;
        public static final int PAIR_UM_PEN_DEFAULT = 17;

        public final int mismatchPenalty;
        public final int gapOpenPenalty;
        public final int gapExtensionPenalty;
        public final int clippingPenalty;
        public final int pairUMPenalty;

        public ScoringScheme(){
            mismatchPenalty = MISMATCH_PEN_DEFAULT;
            gapOpenPenalty = GAP_OPEN_PEN_DEFAULT;
            gapExtensionPenalty = GAP_EXTEN_PEN_DEFAULT;
            clippingPenalty = CLIPPING_PEN_DEFAULT;
            pairUMPenalty = PAIR_UM_PEN_DEFAULT;
        }

        public ScoringScheme(final int mismatchPenalty,
                             final int gapOpenPenalty,
                             final int gapExtensionPenalty,
                             final int clippingPenalty,
                             final int pairUMPenalty){
            this.mismatchPenalty = mismatchPenalty;
            this.gapOpenPenalty = gapOpenPenalty;
            this.gapExtensionPenalty = gapExtensionPenalty;
            this.clippingPenalty = clippingPenalty;
            this.pairUMPenalty = pairUMPenalty;
        }
    }

    /**
     * Preprocess reads before sending them to actual likelihood calculation.
     * Input reads are NOT modified but the returned clones are free to be modified.
     */
    default List<GATKRead> preprocessReads(final List<GATKRead> reads) {
        // TODO: is the deepCopy() call necessary? there's a copy() but doc doesn't seem so assuring...
        return reads.stream().map(GATKRead::deepCopy).collect(Collectors.toList());
    }

    /**
     * Compute the log10 read likelihoods conditional on all possible genotypes in {@code junctions}.
     *
     * @return un-normalized read log10 likelihoods: reads vs all possible alleles, for all samples (see {@link ReadLikelihoods})
     */
    ReadLikelihoods<SVDummyAllele> computeReadLikelihoods(final SampleList sampleList,
                                                          final List<GATKRead> preprocessedReads,
                                                          final Map<String, List<GATKRead>> sample2Reads,
                                                          final SVJunction junctions);
}
