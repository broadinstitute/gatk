package org.broadinstitute.hellbender.tools.picard.analysis.artifacts;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ListMap;
import htsjdk.samtools.util.SequenceUtil;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.SequencingArtifactMetrics.*;

import java.util.*;

/**
 * Keeps track of the AlignmentAccumulators for each artifact / context of interest.
 */
final class ContextAccumulator {

    // are the PE reads expected to face the same direction?
    private final boolean expectedTandemReads;
    private final Map<Transition, Map<String, AlignmentAccumulator>> artifactMap;

    public ContextAccumulator(final Set<String> contexts, final boolean expectedTandemReads) {
        this.expectedTandemReads = expectedTandemReads;
        this.artifactMap = new HashMap<>();
        for (final Transition transition : Transition.values()) {
            this.artifactMap.put(transition, new HashMap<>());
        }
        for (final String context : contexts) {
            final char refBase = getCentralBase(context);
            for (final byte calledBase : SequenceUtil.VALID_BASES_UPPER) {
                final Transition transition = Transition.transitionOf(refBase, (char)calledBase);
                this.artifactMap.get(transition).put(context, new AlignmentAccumulator());
            }
        }
    }

    public void countRecord(final String refContext, final char calledBase, final SAMRecord rec) {
        final char refBase = getCentralBase(refContext);
        final Transition transition = Transition.transitionOf(refBase, calledBase);
        this.artifactMap.get(transition).get(refContext).countRecord(rec);
    }

    /**
     * Core method to compute detailed (i.e. context-by-context) metrics from this accumulator.
     */
    public ListMap<Transition, DetailPair> calculateMetrics(final String sampleAlias, final String library) {
        final ListMap<Transition, DetailPair> detailMetricsMap = new ListMap<>();
        for (final Transition altTransition : Transition.altValues()) {
            final Transition refTransition = altTransition.matchingRef();
            for (final String context : this.artifactMap.get(altTransition).keySet()) {
                // each combination of artifact + context represents a single metric row
                final PreAdapterDetailMetrics preAdapterDetailMetrics = new PreAdapterDetailMetrics();
                final BaitBiasDetailMetrics baitBiasDetailMetrics = new BaitBiasDetailMetrics();

                // populate basic fields
                preAdapterDetailMetrics.SAMPLE_ALIAS = sampleAlias;
                preAdapterDetailMetrics.LIBRARY = library;
                preAdapterDetailMetrics.CONTEXT = context;
                preAdapterDetailMetrics.REF_BASE = altTransition.ref();
                preAdapterDetailMetrics.ALT_BASE = altTransition.call();

                baitBiasDetailMetrics.SAMPLE_ALIAS = sampleAlias;
                baitBiasDetailMetrics.LIBRARY = library;
                baitBiasDetailMetrics.CONTEXT = context;
                baitBiasDetailMetrics.REF_BASE = altTransition.ref();
                baitBiasDetailMetrics.ALT_BASE = altTransition.call();

                // retrieve all the necessary alignment counters.
                final AlignmentAccumulator fwdRefAlignments = this.artifactMap.get(refTransition).get(context);
                final AlignmentAccumulator fwdAltAlignments = this.artifactMap.get(altTransition).get(context);
                final AlignmentAccumulator revRefAlignments = this.artifactMap.get(refTransition.complement()).get(SequenceUtil.reverseComplement(context));
                final AlignmentAccumulator revAltAlignments = this.artifactMap.get(altTransition.complement()).get(SequenceUtil.reverseComplement(context));

                // categorize observations of pre-adapter artifacts
                if (expectedTandemReads) {
                    // if both ends are sequenced on the same strand, then read1/read2 should exhibit the same bias
                    preAdapterDetailMetrics.PRO_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R2_POS + revRefAlignments.R1_NEG + revRefAlignments.R2_NEG;
                    preAdapterDetailMetrics.PRO_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R2_POS + revAltAlignments.R1_NEG + revAltAlignments.R2_NEG;
                    preAdapterDetailMetrics.CON_REF_BASES = fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_NEG + revRefAlignments.R1_POS + revRefAlignments.R2_POS;
                    preAdapterDetailMetrics.CON_ALT_BASES = fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_NEG + revAltAlignments.R1_POS + revAltAlignments.R2_POS;
                } else {
                    // if ends are sequenced on opposite strands, then read1/read2 should exhibit opposite biases
                    preAdapterDetailMetrics.PRO_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R2_NEG + revRefAlignments.R1_NEG + revRefAlignments.R2_POS;
                    preAdapterDetailMetrics.PRO_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R2_NEG + revAltAlignments.R1_NEG + revAltAlignments.R2_POS;
                    preAdapterDetailMetrics.CON_REF_BASES = fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + revRefAlignments.R1_POS + revRefAlignments.R2_NEG;
                    preAdapterDetailMetrics.CON_ALT_BASES = fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + revAltAlignments.R1_POS + revAltAlignments.R2_NEG;
                }

                // categorize observations of bait bias artifacts
                baitBiasDetailMetrics.FWD_CXT_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + fwdRefAlignments.R2_NEG;
                baitBiasDetailMetrics.FWD_CXT_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + fwdAltAlignments.R2_NEG;
                baitBiasDetailMetrics.REV_CXT_REF_BASES = revRefAlignments.R1_POS + revRefAlignments.R1_NEG + revRefAlignments.R2_POS + revRefAlignments.R2_NEG;
                baitBiasDetailMetrics.REV_CXT_ALT_BASES = revAltAlignments.R1_POS + revAltAlignments.R1_NEG + revAltAlignments.R2_POS + revAltAlignments.R2_NEG;

                // calculate error rates + Q-scores
                preAdapterDetailMetrics.calculateDerivedStatistics();
                baitBiasDetailMetrics.calculateDerivedStatistics();

                // add the finalized metrics to the map
                detailMetricsMap.add(altTransition, new DetailPair(preAdapterDetailMetrics, baitBiasDetailMetrics));
            }
        }
        return detailMetricsMap;
    }

    private char getCentralBase(final String context) {
        if (context.length() % 2 == 0) throw new GATKException("Contexts cannot have an even number of bases: " + context);
        else return context.charAt(context.length() / 2);
    }

    /**
     * Little class for breaking down alignments by read1/read2 and positive/negative strand.
     */
    private static class AlignmentAccumulator {
        private long R1_POS = 0;
        private long R1_NEG = 0;
        private long R2_POS = 0;
        private long R2_NEG = 0;

        private void countRecord(final SAMRecord rec) {
            final boolean isNegativeStrand = rec.getReadNegativeStrandFlag();
            final boolean isReadTwo = rec.getReadPairedFlag() && rec.getSecondOfPairFlag();
            if (isReadTwo) {
                if (isNegativeStrand) this.R2_NEG++;
                else this.R2_POS++;
            } else {
                if (isNegativeStrand) this.R1_NEG++;
                else this.R1_POS++;
            }
        }
    }
}
