package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.BreakpointComplications.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.Arrays;
import java.util.List;

/**
 * Based on alignment signature of the input simple chimera, and evidence contig having the chimera, infers
 * <ul>
 *     <li>
 *         exact position of breakpoints following the left-aligning convention,
 *     </li>
 *     <li>
 *         alt haplotype sequence based on given contig sequence
 *     </li>
 *     <li>
 *         complications such as homology, inserted sequence and duplicated ref region, if any.
 *         This inference is actually delegated to {@link BreakpointComplications}.
 *     </li>
 * </ul>
 */
@VisibleForTesting
public abstract class BreakpointsInference {

    protected String upstreamBreakpointRefContig;
    protected String downstreamBreakpointRefContig;
    protected int upstreamBreakpointRefPos;
    protected int downstreamBreakpointRefPos;

    // TODO: 1/26/18 for alt haplotype sequence, we need a single policy to decide if the POS base should be included or not
    protected byte[] altHaplotypeSequence;

    final Tuple2<SimpleInterval, SimpleInterval> getLeftJustifiedBreakpoints() {

        return
                new Tuple2<>(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                             new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos));
    }
    final byte[] getInferredAltHaplotypeSequence() {
        return altHaplotypeSequence;
    }
    abstract BreakpointComplications getComplications();

    protected BreakpointsInference(final SimpleChimera simpleChimera, final byte[] contigSequence, final SAMSequenceDictionary referenceDictionary) {
        resolveComplications(simpleChimera, contigSequence, referenceDictionary);
    }

    static void validateInferredLocations(final SimpleInterval leftBreakpoint,
                                          final SimpleInterval rightBreakpoint,
                                          final SAMSequenceDictionary referenceSequenceDictionary,
                                          final String errorMessage) {

        if ( IntervalUtils.isBefore(rightBreakpoint, leftBreakpoint, referenceSequenceDictionary) ) {
            throw new GATKException.ShouldNeverReachHereException("Inferred novel adjacency reference locations have left location after right location." + errorMessage);
        }

        if ( leftBreakpoint.getEnd() > referenceSequenceDictionary.getSequence(leftBreakpoint.getContig()).getSequenceLength() ) {
            throw new GATKException.ShouldNeverReachHereException("Inferred breakpoint beyond reference sequence length."  + errorMessage);
        }

        if ( rightBreakpoint.getEnd() > referenceSequenceDictionary.getSequence(rightBreakpoint.getContig()).getSequenceLength() ) {
            throw new GATKException.ShouldNeverReachHereException("Inferred breakpoint beyond reference sequence length. " + errorMessage);
        }
    }
    ///////////////
    static BreakpointsInference getInferenceClass(final SimpleChimera simpleChimera,
                                                  final byte[] contigSequence,
                                                  final SAMSequenceDictionary referenceDictionary) {

        switch (simpleChimera.inferType(referenceDictionary)) {
            case INTER_CHR_STRAND_SWITCH_55:
            case INTER_CHR_STRAND_SWITCH_33:
            case INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER:
            case INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER:
                return new InterChromosomeBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
            case INTRA_CHR_REF_ORDER_SWAP:
                return new IntraChrRefOrderSwapBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
            case INTRA_CHR_STRAND_SWITCH_55:
            case INTRA_CHR_STRAND_SWITCH_33:
                if ( simpleChimera.isCandidateInvertedDuplication() ) {
                    return new InvertedDuplicationBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
                } else {
                    return new IntraChrStrandSwitchBreakpointInference(simpleChimera, contigSequence, referenceDictionary);
                }
            case SIMPLE_DEL:
            case RPL:
            case SIMPLE_INS:
                return new SimpleInsertionDeletionBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
            case SMALL_DUP_EXPANSION:
            case DEL_DUP_CONTRACTION:
                return new SmallDuplicationWithPreciseDupRangeBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
            case SMALL_DUP_CPX:
                return new SmallDuplicationWithImpreciseDupRangeBreakpointsInference(simpleChimera, contigSequence, referenceDictionary);
            default:
                throw new GATKException.ShouldNeverReachHereException("Inferred type not recognized for simple chimera:\t" + simpleChimera.toString());
        }
    }

    abstract void resolveComplications(final SimpleChimera simpleChimera, final byte[] contigSequence, final SAMSequenceDictionary referenceDictionary);

    @Override
    public String toString() {
        return "left: " + upstreamBreakpointRefContig + ":" + upstreamBreakpointRefPos +
                " right: " + downstreamBreakpointRefContig + ":" + downstreamBreakpointRefPos +
                " complications: " + getComplications().toString() +
                " alt seq: " + (getInferredAltHaplotypeSequence() == null ? "NULL" : getInferredAltHaplotypeSequence()) ;
    }

    // TODO: 5/28/18 a bug exists here in location inference that when homology exists, see ticket 4883
    //       we are implicitly assuming that the hom. seq. has no micro insertion/deletions
    //       in its alignments at breakpoints (as we are simply adding/subtracting hom. len.)
    //       this is not a serious problem,
    //       unless the net accumulated insertion/deletion length is large
    //       (in the range of tens of bp or longer, since that affects overlap-based evaluation)
    //       similarly, this applies to the duplication cases, where the duplicated copies might not be of exactly the same size
    // =================================================================================================================

    /**
     * Intended to be used for simple insertion, simple deletion, and replacement.
     * NOT FOR SMALL DUPLICATIONS!
     * Alt sequence is
     * <ul>
     *     <li>empty for simple clean deletion,</li>
     *     <li>inserted sequence for inserted sequence</li> todo: see todo above
     *     <li>inserted sequence for replacement</li> todo: see todo above
     * </ul>
     */
    final static class SimpleInsertionDeletionBreakpointsInference extends BreakpointsInference {

        private SimpleInsDelOrReplacementBreakpointComplications complications;

        @Override
        BreakpointComplications getComplications() {
            return complications;
        }

        @Override
        void resolveComplications(final SimpleChimera simpleChimera, final byte[] contigSequence, final SAMSequenceDictionary referenceDictionary) {
            complications =
                    new SimpleInsDelOrReplacementBreakpointComplications(simpleChimera, contigSequence,
                            simpleChimera.firstContigRegionRefSpanAfterSecond(referenceDictionary));
        }


        protected SimpleInsertionDeletionBreakpointsInference (final SimpleChimera simpleChimera,
                                                               final byte[] contigSequence,
                                                               final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);

            upstreamBreakpointRefContig
                    = downstreamBreakpointRefContig
                    = simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig();


            final Tuple2<SimpleInterval, SimpleInterval> coordinateSortedRefSpans = simpleChimera.getCoordinateSortedRefSpans(referenceDictionary);
            final SimpleInterval leftReferenceInterval  = coordinateSortedRefSpans._1,
                                 rightReferenceInterval = coordinateSortedRefSpans._2;
            final int homologyLen = complications.getHomologyForwardStrandRep().length();
            upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
            downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;

            if( complications.insertedSequenceForwardStrandRep.isEmpty() ) { // simple deletion
                altHaplotypeSequence = new byte[0]; // simple deletion has no bases available: (POS, END] is gone!
            } else {
                // note that this will lead to, unless specifically removed, inserted sequence annotation for the deletion call extracted from replacement,
                // even when the inserted sequence is long enough to warrant an insertion,
                // so downstream annotator needs to remove this attribute when linking the INS and DEL
                altHaplotypeSequence = complications.insertedSequenceForwardStrandRep.getBytes();
            }
        }
    }

    ///////////////
    abstract static class SmallDuplicationBreakpointsInference extends BreakpointsInference {

        SmallDuplicationBreakpointsInference(final SimpleChimera simpleChimera,
                                             final byte[] contigSequence,
                                             final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);

            upstreamBreakpointRefContig
                    = downstreamBreakpointRefContig
                    = simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig();
        }
    }

    /**
     * Works for both expansion and contraction.
     * Alt sequence is
     * <ul>
     *     <li>empty if it is contraction</li> todo: see todo above
     *     <li>both copies if it is expansion</li> todo: see todo above
     * </ul>
     */
    final static class SmallDuplicationWithPreciseDupRangeBreakpointsInference extends SmallDuplicationBreakpointsInference {
        private SmallDuplicationWithPreciseDupRangeBreakpointComplications complications;

        @Override
        BreakpointComplications getComplications() {
            return complications;
        }

        @Override
        void resolveComplications(final SimpleChimera simpleChimera, final byte[] contigSequence, final SAMSequenceDictionary referenceDictionary) {
            complications =
                    new SmallDuplicationWithPreciseDupRangeBreakpointComplications(simpleChimera, contigSequence,
                            simpleChimera.firstContigRegionRefSpanAfterSecond(referenceDictionary));
        }

        SmallDuplicationWithPreciseDupRangeBreakpointsInference(final SimpleChimera simpleChimera,
                                                                final byte[] contigSequence,
                                                                final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);

            final int homologyLen = complications.getHomologyForwardStrandRep().length();
            final DistancesBetweenAlignmentsOnRefAndOnRead distances = simpleChimera.getDistancesBetweenAlignmentsOnRefAndOnRead();


            final Tuple2<SimpleInterval, SimpleInterval> coordinateSortedRefSpans = simpleChimera.getCoordinateSortedRefSpans(referenceDictionary);
            final SimpleInterval leftReferenceInterval  = coordinateSortedRefSpans._1,
                                 rightReferenceInterval = coordinateSortedRefSpans._2;

            if ( complications.isDupContraction() ) {
                upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;

                altHaplotypeSequence = new byte[]{};
            } else {
                final int expansionUnitDiff = complications.getDupSeqRepeatNumOnCtg()
                        -
                        complications.getDupSeqRepeatNumOnRef();
                upstreamBreakpointRefPos = leftReferenceInterval.getEnd()
                        - homologyLen
                        - (expansionUnitDiff) * complications.getDupSeqRepeatUnitRefSpan().size();
                downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;

                final List<String> cigarStringsForDupSeqOnCtg = complications.getCigarStringsForDupSeqOnCtgForwardStrandRep();

                final Cigar cigarForFirstCopyOnCtg, cigarForSecondCopyOnCtg;
                if (simpleChimera.isForwardStrandRepresentation) {
                    cigarForFirstCopyOnCtg = TextCigarCodec.decode(cigarStringsForDupSeqOnCtg.get(0));
                    cigarForSecondCopyOnCtg = TextCigarCodec.decode(cigarStringsForDupSeqOnCtg.get(1));
                } else {
                    cigarForFirstCopyOnCtg = TextCigarCodec.decode(cigarStringsForDupSeqOnCtg.get(1));
                    cigarForSecondCopyOnCtg = TextCigarCodec.decode(cigarStringsForDupSeqOnCtg.get(0));
                }

                // TODO: 5/29/18 the small gaps in repeat units will throw this off
                final int zeroBasedStart = distances.firstAlnCtgEnd - cigarForFirstCopyOnCtg.getReadLength();
                final int zeroBasedEnd = distances.secondAlnCtgStart + cigarForSecondCopyOnCtg.getReadLength() - 1;

                altHaplotypeSequence = Arrays.copyOfRange(contigSequence, zeroBasedStart, zeroBasedEnd);
                if (!simpleChimera.isForwardStrandRepresentation) {
                    SequenceUtil.reverseComplement(altHaplotypeSequence);
                }
            }

            validateInferredLocations(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                    new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos),
                    referenceDictionary, toString());
        }
    }

    /**
     * Similar to {@link SmallDuplicationWithPreciseDupRangeBreakpointsInference},
     * except that
     * <ul>
     *     the cigars are not available as the duplicated region is estimated
     *
     * </ul>
     * hence to compensate, alt sequence is
     * <ul>
     *     <li>the left copies if it is contraction</li> todo: see todo above
     *     <li>all copies if it is expansion</li> todo: see todo above
     * </ul>
     */
    final static class SmallDuplicationWithImpreciseDupRangeBreakpointsInference extends SmallDuplicationBreakpointsInference {
        private SmallDuplicationWithImpreciseDupRangeBreakpointComplications complications;

        @Override
        BreakpointComplications getComplications() {
            return complications;
        }


        @Override
        void resolveComplications(final SimpleChimera simpleChimera, final byte[] contigSequence, final SAMSequenceDictionary referenceDictionary) {
            complications =
                    new SmallDuplicationWithImpreciseDupRangeBreakpointComplications(simpleChimera, contigSequence,
                            simpleChimera.firstContigRegionRefSpanAfterSecond(referenceDictionary));
        }

        SmallDuplicationWithImpreciseDupRangeBreakpointsInference(final SimpleChimera simpleChimera,
                                                                  final byte[] contigSequence,
                                                                  final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);

            final DistancesBetweenAlignmentsOnRefAndOnRead distances = simpleChimera.getDistancesBetweenAlignmentsOnRefAndOnRead();

            final int homologyLen = complications.getHomologyForwardStrandRep().length();

            final Tuple2<SimpleInterval, SimpleInterval> coordinateSortedRefSpans = simpleChimera.getCoordinateSortedRefSpans(referenceDictionary);
            final SimpleInterval leftReferenceInterval  = coordinateSortedRefSpans._1,
                                 rightReferenceInterval = coordinateSortedRefSpans._2;
            // PLAN: (note both d_r and d_a < 0; the docs use 1-based coordinate)
            //  For contraction: |d_r| < |d_a|
            //      ALT_SEQ: extract sequence [c2b, c1e], RC if '-'-representation
            //  For expansion:   |d_r| > |d_a|
            //      ALT_SEQ: remember to RC if '-'-representation
            //      ALT_SEQ_BEG: starting at c1e, walk on read backward (with 1st cigar remember to take care of hard clipped bases) a distance corresponding to |d_r| + 2;
            //      ALT_SEQ_END: starting at c2b, walk on read forward  (with 2nd cigar remember to take care of hard clipped bases) a distance corresponding to |d_r| + 2;
            if ( complications.isDupContraction() ) {
                upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;

                altHaplotypeSequence = Arrays.copyOfRange(contigSequence,
                        distances.secondAlnCtgStart - 1, distances.firstAlnCtgEnd);
                if (!simpleChimera.isForwardStrandRepresentation) {
                    SequenceUtil.reverseComplement(altHaplotypeSequence);
                }
            } else {
                final int expansionUnitDiff = complications.getDupSeqRepeatNumOnCtg()
                        -
                        complications.getDupSeqRepeatNumOnRef();
                upstreamBreakpointRefPos = leftReferenceInterval.getEnd()
                        - homologyLen
                        - (expansionUnitDiff) * complications.getDupSeqRepeatUnitRefSpan().size();
                downstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;

                Cigar cigar = simpleChimera.regionWithLowerCoordOnContig.cigarAlong5to3DirectionOfContig;
                int hardClipOffset = cigar.getFirstCigarElement().getOperator().equals(CigarOperator.H) ? cigar.getFirstCigarElement().getLength() : 0;
                int distOnRead = SvCigarUtils
                        .computeAssociatedDistOnRead(cigar,
                                distances.firstAlnCtgEnd - hardClipOffset,
                                -distances.gapBetweenAlignRegionsOnRef,
                                true);
                final int zeroBasedStart = distances.firstAlnCtgEnd - distOnRead;

                cigar = simpleChimera.regionWithHigherCoordOnContig.cigarAlong5to3DirectionOfContig;
                hardClipOffset = cigar.getFirstCigarElement().getOperator().equals(CigarOperator.H) ? cigar.getFirstCigarElement().getLength() : 0;
                distOnRead = SvCigarUtils
                        .computeAssociatedDistOnRead(cigar,
                                distances.secondAlnCtgStart - hardClipOffset,
                                -distances.gapBetweenAlignRegionsOnRef,
                                false);
                final int zeroBasedEnd = distances.secondAlnCtgStart + distOnRead - 1;

                altHaplotypeSequence = Arrays.copyOfRange(contigSequence, zeroBasedStart, zeroBasedEnd);
                if (!simpleChimera.isForwardStrandRepresentation) {
                    SequenceUtil.reverseComplement(altHaplotypeSequence);
                }
            }

            validateInferredLocations(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                    new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos),
                    referenceDictionary, toString());
        }
    }

    ///////////////
    abstract static class BNDTypeBreakpointsInference extends BreakpointsInference {
        protected BNDTypeBreakpointsInference(final SimpleChimera simpleChimera,
                                              final byte[] contigSequence,
                                              final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);
            altHaplotypeSequence = new byte[]{};
        }
    }

    // simple strand-switch breakpoint with alignments NOT overlapping on the read/assembly contig
    final static class IntraChrStrandSwitchBreakpointInference extends BNDTypeBreakpointsInference {

        private IntraChrStrandSwitchBreakpointComplications complications;

        @Override
        BreakpointComplications getComplications() {
            return complications;
        }

        @Override
        void resolveComplications(final SimpleChimera simpleChimera, final byte[] contigSequence, final SAMSequenceDictionary referenceDictionary) {
            complications = new IntraChrStrandSwitchBreakpointComplications(simpleChimera, contigSequence,
                    simpleChimera.firstContigRegionRefSpanAfterSecond(referenceDictionary));
        }

        IntraChrStrandSwitchBreakpointInference(final SimpleChimera simpleChimera,
                                                final byte[] contigSequence,
                                                final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);

            upstreamBreakpointRefContig
                    = downstreamBreakpointRefContig
                    = simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig();

            final int homologyLen = complications.getHomologyForwardStrandRep().length();

            final Tuple2<SimpleInterval, SimpleInterval> coordinateSortedRefSpans = simpleChimera.getCoordinateSortedRefSpans(referenceDictionary);
            final SimpleInterval leftReferenceInterval  = coordinateSortedRefSpans._1,
                                 rightReferenceInterval = coordinateSortedRefSpans._2;
            if (simpleChimera.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE){
                upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                downstreamBreakpointRefPos = rightReferenceInterval.getEnd();
            } else {
                upstreamBreakpointRefPos = leftReferenceInterval.getStart();
                downstreamBreakpointRefPos = rightReferenceInterval.getStart() + homologyLen;
            }

            validateInferredLocations(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                    new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos),
                    referenceDictionary, toString());
        }
    }

    final static class InvertedDuplicationBreakpointsInference extends BNDTypeBreakpointsInference {

        private InvertedDuplicationBreakpointComplications complications;

        @Override
        BreakpointComplications getComplications() {
            return complications;
        }

        @Override
        void resolveComplications(final SimpleChimera simpleChimera, final byte[] contigSequence, final SAMSequenceDictionary referenceDictionary) {
            complications = new InvertedDuplicationBreakpointComplications(simpleChimera, contigSequence,
                    simpleChimera.firstContigRegionRefSpanAfterSecond(referenceDictionary));
        }

        InvertedDuplicationBreakpointsInference(final SimpleChimera simpleChimera,
                                                final byte[] contigSequence,
                                                final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);

            upstreamBreakpointRefContig
                    = downstreamBreakpointRefContig
                    = simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig();

            upstreamBreakpointRefPos = complications.getDupSeqRepeatUnitRefSpan().getStart() - 1;
            downstreamBreakpointRefPos = complications.getDupSeqRepeatUnitRefSpan().getEnd();

            altHaplotypeSequence = extractAltHaplotypeForInvDup(simpleChimera, contigSequence);

            validateInferredLocations(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                    new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos),
                    referenceDictionary, toString());
        }

        private static byte[] extractAltHaplotypeForInvDup(final SimpleChimera simpleChimera, final byte[] contigSeq) {

            final AlignmentInterval firstAlignmentInterval  = simpleChimera.regionWithLowerCoordOnContig;
            final AlignmentInterval secondAlignmentInterval = simpleChimera.regionWithHigherCoordOnContig;

            final int start, end; // intended to be 0-based, semi-open [start, end)
            final boolean needRC;
            // below we need to use cigars of the provided alignments to compute how long we need to walk on the read
            // so that we can "start" to or "end" to collect bases for alternative haplotype sequence,
            // because one could imagine either alignment has long flanking region that is far from the affected reference region.
            if (firstAlignmentInterval.forwardStrand) {
                final int alpha = firstAlignmentInterval.referenceSpan.getStart(),
                        omega = secondAlignmentInterval.referenceSpan.getStart();
                if (alpha <= omega) {
                    final int walkOnReadUntilDuplicatedSequence ;
                    if (alpha == omega) {
                        walkOnReadUntilDuplicatedSequence = 0;
                    } else {
                        walkOnReadUntilDuplicatedSequence = SvCigarUtils.computeAssociatedDistOnRead(firstAlignmentInterval.cigarAlong5to3DirectionOfContig,
                                firstAlignmentInterval.startInAssembledContig, omega - alpha, false);
                    }
                    start = firstAlignmentInterval.startInAssembledContig + walkOnReadUntilDuplicatedSequence - 1;
                    end = secondAlignmentInterval.endInAssembledContig;
                    needRC = false;
                } else {
                    final int walkOnReadUntilDuplicatedSequence = SvCigarUtils.computeAssociatedDistOnRead(secondAlignmentInterval.cigarAlong5to3DirectionOfContig,
                            secondAlignmentInterval.endInAssembledContig, alpha - omega, true);
                    start = firstAlignmentInterval.startInAssembledContig - 1;
                    end = secondAlignmentInterval.endInAssembledContig - walkOnReadUntilDuplicatedSequence;
                    needRC = true;
                }
            } else {
                final int alpha = firstAlignmentInterval.referenceSpan.getEnd(),
                        omega = secondAlignmentInterval.referenceSpan.getEnd();
                if (alpha >= omega) {
                    final int walkOnReadUntilDuplicatedSequence ;
                    if (alpha == omega) {
                        walkOnReadUntilDuplicatedSequence = 0;
                    } else {
                        walkOnReadUntilDuplicatedSequence = SvCigarUtils.computeAssociatedDistOnRead(firstAlignmentInterval.cigarAlong5to3DirectionOfContig,
                                firstAlignmentInterval.startInAssembledContig, alpha - omega, false);
                    }
                    start = firstAlignmentInterval.startInAssembledContig + walkOnReadUntilDuplicatedSequence - 1;
                    end = secondAlignmentInterval.endInAssembledContig;
                    needRC = true;
                } else {
                    final int walkOnReadUntilDuplicatedSequence = SvCigarUtils.computeAssociatedDistOnRead(secondAlignmentInterval.cigarAlong5to3DirectionOfContig,
                            secondAlignmentInterval.endInAssembledContig, omega - alpha, true);
                    start = firstAlignmentInterval.startInAssembledContig - 1;
                    end = secondAlignmentInterval.endInAssembledContig - walkOnReadUntilDuplicatedSequence;
                    needRC = false;
                }
            }

            final byte[] seq = Arrays.copyOfRange(contigSeq, start, end);
            if (needRC) SequenceUtil.reverseComplement(seq, 0, seq.length);
            return seq;
        }
    }

    final static class IntraChrRefOrderSwapBreakpointsInference extends BNDTypeBreakpointsInference {

        private IntraChrRefOrderSwapBreakpointComplications complications;

        @Override
        BreakpointComplications getComplications() {
            return complications;
        }

        @Override
        void resolveComplications(final SimpleChimera simpleChimera, final byte[] contigSequence, final SAMSequenceDictionary referenceDictionary) {
            complications = new IntraChrRefOrderSwapBreakpointComplications(simpleChimera, contigSequence,
                    simpleChimera.firstContigRegionRefSpanAfterSecond(referenceDictionary));
        }

        IntraChrRefOrderSwapBreakpointsInference(final SimpleChimera simpleChimera,
                                                 final byte[] contigSequence,
                                                 final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);

            final int homologyLen = complications.getHomologyForwardStrandRep().length();
            final Tuple2<SimpleInterval, SimpleInterval> coordinateSortedRefSpans = simpleChimera.getCoordinateSortedRefSpans(referenceDictionary);
            final SimpleInterval leftRefSpan  = coordinateSortedRefSpans._1,
                                 rightRefSpan = coordinateSortedRefSpans._2;
            upstreamBreakpointRefContig
                    = downstreamBreakpointRefContig
                    = simpleChimera.regionWithLowerCoordOnContig.referenceSpan.getContig();
            upstreamBreakpointRefPos = leftRefSpan.getStart();
            downstreamBreakpointRefPos = rightRefSpan.getEnd() - homologyLen;

            validateInferredLocations(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                    new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos),
                    referenceDictionary, toString());
        }
    }

    /**
     * For computing exact and left-adjusted breakpoint locations of inter-chromosome novel adjacency,
     * with or without strand switch.
     */
    final static class InterChromosomeBreakpointsInference extends BNDTypeBreakpointsInference {

        private InterChromosomeBreakpointComplications complications;

        @Override
        BreakpointComplications getComplications() {
            return complications;
        }

        @Override
        void resolveComplications(final SimpleChimera simpleChimera, final byte[] contigSequence,
                                  final SAMSequenceDictionary referenceDictionary) {
            complications = new InterChromosomeBreakpointComplications(simpleChimera, contigSequence,
                    simpleChimera.firstContigRegionRefSpanAfterSecond(referenceDictionary));
        }

        InterChromosomeBreakpointsInference(final SimpleChimera simpleChimera,
                                            final byte[] contigSequence,
                                            final SAMSequenceDictionary referenceDictionary) {
            super(simpleChimera, contigSequence, referenceDictionary);

            determineRefContigs(simpleChimera, referenceDictionary);

            extractRefPositions(simpleChimera, complications, referenceDictionary);

            validateInferredLocations(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                    new SimpleInterval(downstreamBreakpointRefContig, downstreamBreakpointRefPos, downstreamBreakpointRefPos),
                    referenceDictionary, toString());
        }

        private void extractRefPositions(final SimpleChimera ca,
                                         final InterChromosomeBreakpointComplications complication,
                                         final SAMSequenceDictionary referenceDictionary) {
            final int homologyLen = complication.getHomologyForwardStrandRep().length();
            final boolean firstInPartner = isFirstInPartner(ca, referenceDictionary);
            if (firstInPartner) {
                switch (ca.strandSwitch) {
                    case NO_SWITCH:
                        if (ca.isForwardStrandRepresentation) {
                            upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            downstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                        } else {
                            upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                            downstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                        }
                        break;
                    case FORWARD_TO_REVERSE:
                        upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                        downstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd();
                        break;
                    case REVERSE_TO_FORWARD:
                        upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                        downstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart() + homologyLen;
                        break;
                    default: throw new GATKException.ShouldNeverReachHereException("Unseen strand switch case for: " + ca.toString());
                }
            } else {
                switch (ca.strandSwitch) {
                    case NO_SWITCH:
                        if (ca.isForwardStrandRepresentation) {
                            upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                            downstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                        } else {
                            upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            downstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                        }
                        break;
                    case FORWARD_TO_REVERSE:
                        upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                        downstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd();
                        break;
                    case REVERSE_TO_FORWARD:
                        upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                        downstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart() + homologyLen;
                        break;
                    default: throw new GATKException.ShouldNeverReachHereException("Unseen strand switch case for: " + ca.toString());
                }
            }
        }

        private void determineRefContigs(SimpleChimera ca, SAMSequenceDictionary referenceDictionary) {
            final boolean firstInPartner = isFirstInPartner(ca, referenceDictionary);
            if (firstInPartner) {
                upstreamBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
                downstreamBreakpointRefContig = ca.regionWithHigherCoordOnContig.referenceSpan.getContig();
            } else {
                upstreamBreakpointRefContig = ca.regionWithHigherCoordOnContig.referenceSpan.getContig();
                downstreamBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
            }
        }

        private static boolean isFirstInPartner(final SimpleChimera ca, final SAMSequenceDictionary referenceDictionary) {
            switch (ca.strandSwitch) {
                case NO_SWITCH: return 0 > IntervalUtils.compareContigs(ca.regionWithLowerCoordOnContig.referenceSpan,
                        ca.regionWithHigherCoordOnContig.referenceSpan, referenceDictionary);
                case FORWARD_TO_REVERSE: case REVERSE_TO_FORWARD:
                    return ca.isForwardStrandRepresentation;
                default:
                    throw new GATKException.ShouldNeverReachHereException("Unseen strand switch case for: " + ca.toString());
            }
        }
    }
}
