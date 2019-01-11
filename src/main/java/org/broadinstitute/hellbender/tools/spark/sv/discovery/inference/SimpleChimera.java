package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.List;

/**
 * Conceptually, a simple chimera represents the junction on
 * {@link AssemblyContigWithFineTunedAlignments} that have
 * exactly two good alignments.
 */
@DefaultSerializer(SimpleChimera.Serializer.class)
public class SimpleChimera {

    public final String sourceContigName;

    public final AlignmentInterval regionWithLowerCoordOnContig;
    public final AlignmentInterval regionWithHigherCoordOnContig;

    final StrandSwitch strandSwitch;
    final boolean isForwardStrandRepresentation;

    public final List<String> insertionMappings;

    public final String goodNonCanonicalMappingSATag;

    // =================================================================================================================
    // (conditional) construction block

    protected SimpleChimera(final Kryo kryo, final Input input) {

        this.sourceContigName = input.readString();

        this.regionWithLowerCoordOnContig = kryo.readObject(input, AlignmentInterval.class);
        this.regionWithHigherCoordOnContig = kryo.readObject(input, AlignmentInterval.class);

        this.strandSwitch = StrandSwitch.values()[input.readInt()];
        this.isForwardStrandRepresentation = input.readBoolean();

        final int insertionsMappingSize = input.readInt();
        this.insertionMappings = new ArrayList<>( insertionsMappingSize );
        for(int i = 0; i < insertionsMappingSize; ++i) {
            insertionMappings.add(input.readString());
        }
        goodNonCanonicalMappingSATag = input.readString();
    }

    @VisibleForTesting
    public SimpleChimera(final String sourceContigName,
                         final AlignmentInterval regionWithLowerCoordOnContig, final AlignmentInterval regionWithHigherCoordOnContig,
                         final StrandSwitch strandSwitch, final boolean isForwardStrandRepresentation,
                         final List<String> insertionMappings, final String goodNonCanonicalMappingSATag) {
        this.sourceContigName = sourceContigName;
        this.regionWithLowerCoordOnContig = regionWithLowerCoordOnContig;
        this.regionWithHigherCoordOnContig = regionWithHigherCoordOnContig;
        this.strandSwitch = strandSwitch;
        this.isForwardStrandRepresentation = isForwardStrandRepresentation;
        this.insertionMappings = insertionMappings;
        this.goodNonCanonicalMappingSATag = goodNonCanonicalMappingSATag;
    }

    /**
     * Construct a new SimpleChimera from two alignment intervals.
     * Assumes {@code intervalWithLowerCoordOnContig} has a lower {@link AlignmentInterval#startInAssembledContig}
     * than {@code regionWithHigherCoordOnContig}.
     */
    public SimpleChimera(final AlignmentInterval intervalWithLowerCoordOnContig, final AlignmentInterval intervalWithHigherCoordOnContig,
                         final List<String> insertionMappings, final String sourceContigName, final String goodNonCanonicalMappingSATag,
                         final SAMSequenceDictionary referenceDictionary) {

        this.sourceContigName = sourceContigName;

        this.regionWithLowerCoordOnContig = intervalWithLowerCoordOnContig;
        this.regionWithHigherCoordOnContig = intervalWithHigherCoordOnContig;

        this.strandSwitch = determineStrandSwitch(intervalWithLowerCoordOnContig, intervalWithHigherCoordOnContig);

        this.isForwardStrandRepresentation =
                isForwardStrandRepresentation(intervalWithLowerCoordOnContig, intervalWithHigherCoordOnContig, strandSwitch, referenceDictionary);

        this.insertionMappings = insertionMappings;
        this.goodNonCanonicalMappingSATag = goodNonCanonicalMappingSATag;
    }

    /**
     * Roughly similar to
     * DiscoverVariantsFromContigAlignmentsSAMSpark#nextAlignmentMayBeInsertion(AlignmentInterval, AlignmentInterval, Integer, Integer, boolean):
     *  1) either alignment may have very low mapping quality (a more relaxed mapping quality threshold);
     *  2) either alignment may consume only a "short" part of the contig, or if assuming that the alignment consumes
     *     roughly the same amount of ref bases and read bases, has isAlignment that is too short
     */
    static boolean splitPairStrongEnoughEvidenceForCA(final AlignmentInterval intervalOne,
                                                      final AlignmentInterval intervalTwo,
                                                      final int mapQThresholdInclusive,
                                                      final int alignmentLengthThresholdInclusive) {

        if (intervalOne.mapQual < mapQThresholdInclusive || intervalTwo.mapQual < mapQThresholdInclusive)
            return false;

        // TODO: 2/2/18 improve annotation for alignment length: compared to #firstAlignmentIsTooShort(),
        // we are not subtracting alignments' overlap on the read, i.e. we are not filtering alignments based on their unique read span size,
        // but downstream analysis should have this information via an annotation, the current annotation is not up for this task
        return Math.min(intervalOne.getSizeOnRead(), intervalTwo.getSizeOnRead()) >= alignmentLengthThresholdInclusive;
    }

    /**
     * An SV event could be detected from a contig that seem to originate from the forward or reverse strand of the reference,
     * besides the annotation that the alignment flanking regions might flank either side of the two breakpoints.
     *
     * <p> The definition for '+' representation is such that:
     *     <ul>
     *         <li>For events involving alignments to the same chromosome (with or without strand switch),
     *              1) the two alignments are mapped to the '+' reference strand when there is NO strand switch
     *              2) {@code regionWithLowerCoordOnContig} ends earlier on the reference than
     *                 {@code regionWithHigherCoordOnContig} does when
     *                 {@code regionWithLowerCoordOnContig} is '+' and {@code regionWithHigherCoordOnContig} is '-'
     *              3) {@code regionWithLowerCoordOnContig} starts earlier on the reference than
     *                 {@code regionWithHigherCoordOnContig} does when
     *                 {@code regionWithLowerCoordOnContig} is '-' and {@code regionWithHigherCoordOnContig} is '+'
     *         </li>
     *         <li>For events involving alignments to different chromosomes WITHOUT strand switch,
     *              the two alignments are mapped to the '+' reference strand
     *         </li>
     *         <li>For events involving alignments to different chromosomes with strand switch,
     *              {@code regionWithLowerCoordOnContig} is a mapping to a chromosome with lower index, according to
     *              {@code referenceDictionary}
     *         </li>
     *     </ul>
     * </p>
     *
     */
    @VisibleForTesting
    static boolean isForwardStrandRepresentation(final AlignmentInterval regionWithLowerCoordOnContig,
                                                 final AlignmentInterval regionWithHigherCoordOnContig,
                                                 final StrandSwitch strandSwitch,
                                                 final SAMSequenceDictionary referenceDictionary) {

        final boolean mappedToSameChr = regionWithLowerCoordOnContig.referenceSpan.getContig()
                .equals(regionWithHigherCoordOnContig.referenceSpan.getContig());
        if (mappedToSameChr) {
            switch (strandSwitch) {
                case NO_SWITCH: return regionWithLowerCoordOnContig.forwardStrand;
                case FORWARD_TO_REVERSE: return regionWithLowerCoordOnContig.referenceSpan.getEnd() < regionWithHigherCoordOnContig.referenceSpan.getEnd();
                case REVERSE_TO_FORWARD: return regionWithLowerCoordOnContig.referenceSpan.getStart() < regionWithHigherCoordOnContig.referenceSpan.getStart();
                default: throw new IllegalArgumentException("Seeing unexpected strand switch case: " + strandSwitch.name());
            }
        } else {
            if (strandSwitch == StrandSwitch.NO_SWITCH) {
                return regionWithLowerCoordOnContig.forwardStrand;
            } else {
                return IntervalUtils.compareContigs(regionWithLowerCoordOnContig.referenceSpan,
                        regionWithHigherCoordOnContig.referenceSpan, referenceDictionary)
                        < 0;
            }
        }
    }

    @VisibleForTesting
    static StrandSwitch determineStrandSwitch(final AlignmentInterval first, final AlignmentInterval second) {
        if (first.forwardStrand == second.forwardStrand) {
            return StrandSwitch.NO_SWITCH;
        } else {
            return first.forwardStrand ? StrandSwitch.FORWARD_TO_REVERSE : StrandSwitch.REVERSE_TO_FORWARD;
        }
    }


    /**
     * For the two alignments, test if the alignment that has lower coordinate on read/contig
     * has a higher coordinate (as specified by {@code referenceDictionary}) on reference.
     *
     * This could happen for simple insertions/deletions when the evidence contig is
     * <ul>
     *     <li>a reverse strand representation of a simple event, or</li>
     *     <li>a reverse strand representation of left breakpoint of inversion, or</li>
     *     <li>a reverse strand representation of right breakpoint of inversion, or</li>
     *     <li>a forward strand representation of a ref-order swap event (e.g. large tandem duplication breakpoint suspect)</li>
     *     <li>a reverse strand representation of inter-chromosomal events, except one scenario
     *          where the two forward strand mapping alignments jump from a higher ref coordinate to a lower coordinate
     *          (see definition/convention of strand representation in
     *          {@link #isForwardStrandRepresentation(AlignmentInterval, AlignmentInterval, StrandSwitch, SAMSequenceDictionary)}
     *     </li>
     * </ul>
     */
    boolean firstContigRegionRefSpanAfterSecond(final SAMSequenceDictionary referenceDictionary){
        return IntervalUtils.compareLocatables(regionWithLowerCoordOnContig.referenceSpan, regionWithHigherCoordOnContig.referenceSpan,
                                               referenceDictionary) > 0;
    }

    Tuple2<SimpleInterval, SimpleInterval> getCoordinateSortedRefSpans(final SAMSequenceDictionary referenceDictionary) {

        if (firstContigRegionRefSpanAfterSecond(referenceDictionary)) {
            return new Tuple2<>(regionWithHigherCoordOnContig.referenceSpan, regionWithLowerCoordOnContig.referenceSpan);
        } else {
            return new Tuple2<>(regionWithLowerCoordOnContig.referenceSpan, regionWithHigherCoordOnContig.referenceSpan);
        }
    }

    // =================================================================================================================

    /**
     * Implementing a logic, where based on the simple chimera, which show's how a read or assembly contig's
     * alignments overlap or are distant from each other, infer the possible simple breakpoint type.
     * @throws IllegalArgumentException when the simple chimera indicates strand switch or simple translocation or incomplete picture.
     */
    TypeInferredFromSimpleChimera inferType(final SAMSequenceDictionary referenceDictionary) {

        if ( isCandidateSimpleTranslocation()) { // see {@link SimpleChimera.isCandidateSimpleTranslocation()} for definition
            final boolean sameChromosomeEvent =
                    regionWithLowerCoordOnContig.referenceSpan.getContig()
                            .equals(regionWithHigherCoordOnContig.referenceSpan.getContig());
            if ( sameChromosomeEvent ) {
                return TypeInferredFromSimpleChimera.INTRA_CHR_REF_ORDER_SWAP;
            } else {
                if (strandSwitch.equals(StrandSwitch.FORWARD_TO_REVERSE)) {
                    return TypeInferredFromSimpleChimera.INTER_CHR_STRAND_SWITCH_55;
                } else if (strandSwitch.equals(StrandSwitch.REVERSE_TO_FORWARD)) {
                    return TypeInferredFromSimpleChimera.INTER_CHR_STRAND_SWITCH_33;
                } else { // no switch, but still need to distinguish between cases of pair WX vs UV in Fig. 7 in Section 5.4 of VCF spec ver.4.2
                    if (isForwardStrandRepresentation != firstContigRegionRefSpanAfterSecond(referenceDictionary) ){
                        return TypeInferredFromSimpleChimera.INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER;
                    } else {
                        return TypeInferredFromSimpleChimera.INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER;
                    }
                }
            }
        } else {
            if (strandSwitch.equals(StrandSwitch.FORWARD_TO_REVERSE)) {
                // TODO: 9/9/17 the case involves an inversion, could be retired once same chr strand-switch BND calls are evaluated.
                return TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_55;
            } else if (strandSwitch.equals(StrandSwitch.REVERSE_TO_FORWARD)) {
                return TypeInferredFromSimpleChimera.INTRA_CHR_STRAND_SWITCH_33;
            } else {
                final DistancesBetweenAlignmentsOnRefAndOnRead distances = getDistancesBetweenAlignmentsOnRefAndOnRead();
                final int distBetweenAlignRegionsOnRef = distances.gapBetweenAlignRegionsOnRef, // distance-1 between the two regions on reference, denoted as d1 in the comments below
                          distBetweenAlignRegionsOnCtg = distances.gapBetweenAlignRegionsOnCtg; // distance-1 between the two regions on contig, denoted as d2 in the comments below
                if (distBetweenAlignRegionsOnRef > 0) {        // Deletion:
                    if (distBetweenAlignRegionsOnCtg <= 0) {     // simple deletion when == 0; with homology when < 0
                        return TypeInferredFromSimpleChimera.SIMPLE_DEL;
                    } else {
                        return TypeInferredFromSimpleChimera.RPL;
                    }
                } else if (distBetweenAlignRegionsOnRef < 0) {
                    if (distBetweenAlignRegionsOnCtg >= 0) { // Tandem repeat expansion:   reference bases [r1e-|d1|+1, r1e] to contig bases [c1e-|d1|+1, c1e] and [c2b, c2b+|d1|-1] with optional inserted sequence [c1e+1, c2b-1] in between the two intervals on contig
                        return TypeInferredFromSimpleChimera.SMALL_DUP_EXPANSION;
                    } else {  // complicated case, see below
                        // Deletion:  duplication with repeat number N1 on reference, N2 on contig, such that N1 <= 2*N2 (and N2<N1);
                        // Insertion: duplication with repeat number N1 on reference, N2 on contig, such that N2 <= 2*N1 (and N1<N2);
                        // in both cases, the equal sign on the right can be taken only when there's pseudo-homology between starting bases of the duplicated sequence and starting bases of the right flanking region
                        // the reference system with a shorter overlap (i.e. with less-negative distance between regions) has a higher repeat number
                        return TypeInferredFromSimpleChimera.SMALL_DUP_CPX;
                    }
                } else {  // gapBetweenAlignRegionsOnRef == 0
                    if (distBetweenAlignRegionsOnCtg > 0) { // Insertion: simple insertion, inserted sequence is the sequence [c1e+1, c2b-1] on the contig
                        return TypeInferredFromSimpleChimera.SIMPLE_INS;
                    } else if (distBetweenAlignRegionsOnCtg < 0) { // Tandem repeat contraction: reference has two copies but one copy was deleted on the contig; duplicated sequence on reference are [r1e-|d2|+1, r1e] and [r2b, r2b+|d2|-1]
                        return TypeInferredFromSimpleChimera.DEL_DUP_CONTRACTION;
                    } else { // both == 0 => SNP & indel
                        throw new GATKException.ShouldNeverReachHereException(
                                "Detected badly parsed chimeric alignment for identifying SV breakpoints; no rearrangement found: "
                                        + toString());
                    }
                }
            }
        }
    }

    public boolean isNeitherIncompleteNorSimpleTranslocation() {
        if ( hasIncompletePicture() )
            return false;
        return ! isCandidateSimpleTranslocation();
    }

    /**
     * Determine if the chimeric alignment indicates a simple translocation.
     * Simple translocations are defined here and at this time as:
     * <ul>
     *     <li>
     *         evidence alignments doesn't show signature specified in {@link #hasIncompletePicture()}
     *     </li>
     *     <li>inter-chromosomal translocations, i.e. novel adjacency between different reference chromosomes, or</li>
     *      <li>intra-chromosomal translocation that DOES NOT involve a strand switch, that is:
     *          novel adjacency between reference locations on the same chromosome involving reference-order switch but NO strand switch,
     *          but in the meantime, the two inducing alignments CANNOT overlap each other since that would point to
     *          incomplete picture, hence not "simple" anymore.
     *      </li>
     * </ul>
     *
     * <p>
     *     Note that the above definition specifically does not cover the case where
     *     the suggested novel adjacency are linking two reference locations that are
     *     on the same chromosome involving a strand switch.
     *     The reason is: to resolve/interpret such cases
     *     (distinguish between a real inversion or dispersed inverted duplication),
     *     we need other types of evidence in addition to contig alignment signatures.
     * </p>
     */
    @VisibleForTesting
    boolean isCandidateSimpleTranslocation() {

        if ( hasIncompletePicture() ) // note that this implies inter-chromosome events are kept
            return false;

        if (!regionWithLowerCoordOnContig.referenceSpan.getContig()
                .equals(regionWithHigherCoordOnContig.referenceSpan.getContig()))
            return true;

        if ( !strandSwitch.equals(StrandSwitch.NO_SWITCH) ) // note that simple chimera having complete picture could have strand switch
            return false;

        final SimpleInterval referenceSpanOne = regionWithLowerCoordOnContig.referenceSpan,
                             referenceSpanTwo = regionWithHigherCoordOnContig.referenceSpan;

        if (regionWithLowerCoordOnContig.forwardStrand) {
            return referenceSpanOne.getStart() > referenceSpanTwo.getEnd();
        } else {
            return referenceSpanTwo.getStart() > referenceSpanOne.getEnd();
        }
    }

    /**
     * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromTwoAlignments(AlignmentInterval, AlignmentInterval)}.
     */
    private boolean hasIncompletePicture() {
        return AssemblyContigWithFineTunedAlignments.hasIncompletePictureFromTwoAlignments(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig);
    }

    // TODO: 5/5/18 Note that the use of the following predicate is currently obsoleted. Detail below: by
    //      This predicate is currently used in two places (excluding appearance in comments):
    //      `BreakpointComplications.IntraChrStrandSwitchBreakpointComplications`,
    //      where it is use to test if the input SimpleChimera indicates an inverse tandem duplication and trigger the logic for inferring duplicated region; and
    //      `BreakpointsInference.IntraChrStrandSwitchBreakpointInference`, where it is used for breakpoints inference.
    //      The problem is, the contig holding this simple chimera will NOT even be sent here,
    //      because {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromTwoAlignments()} defines
    //      a contig with simple chimera that has strand switch and the two alignments overlaps on reference
    //      as "incomplete" (the duplicated region is not complete),
    //      so in practice the two uses are NOT going to be triggered.
    //      But when we come back later and see what can be extracted from such "incomplete" contigs,
    //      these code could be useful again. So it is kept.
    // another old todo: see ticket #3529 (Change to a more principled criterion than more than half of alignments overlapping)
    /**
     * @return true iff the two alignments of the assembly contig are
     *         1) mappings to different strands on the same chromosome, and
     *         2) overlapping on reference is more than half of the two AI's minimal read span.
     */
    @VisibleForTesting
    boolean isCandidateInvertedDuplication() {
        if (regionWithLowerCoordOnContig.forwardStrand == regionWithHigherCoordOnContig.forwardStrand)
            return false;
        return 2 * AlignmentInterval.overlapOnRefSpan(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig)
                >
                Math.min(regionWithLowerCoordOnContig.getSizeOnRead(), regionWithHigherCoordOnContig.getSizeOnRead());
    }

    /**
     * Struct to represent the (distance - 1) between boundaries of the two alignments represented by this CA,
     * on reference, and on read.
     * For example,
     * two alignments have ref spans  1:100-200, 1:151-250
     *                     read spans 1:100, 81-181
     * then their distance on reference would be -50, and on read would be -20.
     *
     * Note that
     *
     * Note that this concept is ONLY applicable to chimeric alignments that are
     * {@link #isNeitherIncompleteNorSimpleTranslocation()} and
     * {@link #determineStrandSwitch(AlignmentInterval, AlignmentInterval)} == {@link StrandSwitch#NO_SWITCH}
     */
    public static final class DistancesBetweenAlignmentsOnRefAndOnRead {
        final int gapBetweenAlignRegionsOnRef; // size of the gap between the two regions on reference, could be negative
        final int gapBetweenAlignRegionsOnCtg; // size of the gap between the two regions on contig, could be negative

        final int leftAlnRefEnd;
        final int rightAlnRefStart;
        final int firstAlnCtgEnd;
        final int secondAlnCtgStart;

        public DistancesBetweenAlignmentsOnRefAndOnRead(final int gapBetweenAlignRegionsOnRef,
                                                        final int gapBetweenAlignRegionsOnCtg,
                                                        final int leftAlnRefEnd,
                                                        final int rightAlnRefStart,
                                                        final int firstAlnCtgEnd,
                                                        final int secondAlnCtgStart) {
            this.gapBetweenAlignRegionsOnRef = gapBetweenAlignRegionsOnRef;
            this.gapBetweenAlignRegionsOnCtg = gapBetweenAlignRegionsOnCtg;
            this.leftAlnRefEnd = leftAlnRefEnd;
            this.rightAlnRefStart = rightAlnRefStart;
            this.firstAlnCtgEnd = firstAlnCtgEnd;
            this.secondAlnCtgStart = secondAlnCtgStart;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final DistancesBetweenAlignmentsOnRefAndOnRead that = (DistancesBetweenAlignmentsOnRefAndOnRead) o;

            if (gapBetweenAlignRegionsOnRef != that.gapBetweenAlignRegionsOnRef) return false;
            if (gapBetweenAlignRegionsOnCtg != that.gapBetweenAlignRegionsOnCtg) return false;
            if (leftAlnRefEnd != that.leftAlnRefEnd) return false;
            if (rightAlnRefStart != that.rightAlnRefStart) return false;
            if (firstAlnCtgEnd != that.firstAlnCtgEnd) return false;
            return secondAlnCtgStart == that.secondAlnCtgStart;
        }

        @Override
        public int hashCode() {
            int result = gapBetweenAlignRegionsOnRef;
            result = 31 * result + gapBetweenAlignRegionsOnCtg;
            result = 31 * result + leftAlnRefEnd;
            result = 31 * result + rightAlnRefStart;
            result = 31 * result + firstAlnCtgEnd;
            result = 31 * result + secondAlnCtgStart;
            return result;
        }

        @Override
        public String toString() {
            final StringBuilder sb = new StringBuilder("DistancesBetweenAlignmentsOnRefAndOnRead{");
            sb.append("gapBetweenAlignRegionsOnRef=").append(gapBetweenAlignRegionsOnRef);
            sb.append(", gapBetweenAlignRegionsOnCtg=").append(gapBetweenAlignRegionsOnCtg);
            sb.append(", leftAlnRefEnd=").append(leftAlnRefEnd);
            sb.append(", rightAlnRefStart=").append(rightAlnRefStart);
            sb.append(", firstAlnCtgEnd=").append(firstAlnCtgEnd);
            sb.append(", secondAlnCtgStart=").append(secondAlnCtgStart);
            sb.append('}');
            return sb.toString();
        }
    }

    /**
     * @throws IllegalArgumentException when the simple chimera is either incomplete or simple translocation, or
     *                                  when the chimera involves strand switch
     */
    DistancesBetweenAlignmentsOnRefAndOnRead getDistancesBetweenAlignmentsOnRefAndOnRead() {
        if ( ! (isNeitherIncompleteNorSimpleTranslocation() && strandSwitch.equals(StrandSwitch.NO_SWITCH)) ) {
            throw new UnsupportedOperationException(
                    "Assumption that the simple chimera is neither incomplete picture nor simple translocation is violated.\n" +
                            toString());
        }
        final AlignmentInterval firstContigRegion  = regionWithLowerCoordOnContig;
        final AlignmentInterval secondContigRegion = regionWithHigherCoordOnContig;
        final SimpleInterval leftReferenceSpan, rightReferenceSpan;
        if (isForwardStrandRepresentation) {
            leftReferenceSpan = firstContigRegion.referenceSpan;
            rightReferenceSpan = secondContigRegion.referenceSpan;
        } else {
            leftReferenceSpan = secondContigRegion.referenceSpan;
            rightReferenceSpan = firstContigRegion.referenceSpan;
        }

        final int r1e = leftReferenceSpan.getEnd(),
                r2b = rightReferenceSpan.getStart(),
                c1e = firstContigRegion.endInAssembledContig,
                c2b = secondContigRegion.startInAssembledContig;

        return new DistancesBetweenAlignmentsOnRefAndOnRead(r2b - r1e - 1,
                c2b - c1e - 1,
                r1e,
                r2b,
                c1e,
                c2b);
    }

    // =================================================================================================================

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("SimpleChimera{");
        sb.append("sourceContigName='").append(sourceContigName).append('\'');
        sb.append(", regionWithLowerCoordOnContig=").append(regionWithLowerCoordOnContig);
        sb.append(", regionWithHigherCoordOnContig=").append(regionWithHigherCoordOnContig);
        sb.append(", strandSwitch=").append(strandSwitch);
        sb.append(", isForwardStrandRepresentation=").append(isForwardStrandRepresentation);
        sb.append(", insertionMappings=").append(insertionMappings);
        sb.append(", goodNonCanonicalMappingSATag='").append(goodNonCanonicalMappingSATag).append('\'');
        sb.append('}');
        return sb.toString();
    }

    protected void serialize(final Kryo kryo, final Output output) {

        output.writeString(sourceContigName);

        kryo.writeObject(output, regionWithLowerCoordOnContig);
        kryo.writeObject(output, regionWithHigherCoordOnContig);

        output.writeInt(strandSwitch.ordinal());
        output.writeBoolean(isForwardStrandRepresentation);

        final int insertionsMappingSize = insertionMappings.size();
        output.writeInt(insertionsMappingSize);
        insertionMappings.forEach(output::writeString);
        output.writeString(goodNonCanonicalMappingSATag);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SimpleChimera> {
        @Override
        public void write(final Kryo kryo, final Output output, final SimpleChimera simpleChimera) {
            simpleChimera.serialize(kryo, output);
        }

        @Override
        public SimpleChimera read(final Kryo kryo, final Input input, final Class<SimpleChimera> klass ) {
            return new SimpleChimera(kryo, input);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SimpleChimera that = (SimpleChimera) o;

        if (isForwardStrandRepresentation != that.isForwardStrandRepresentation) return false;
        if (!sourceContigName.equals(that.sourceContigName)) return false;
        if (!regionWithLowerCoordOnContig.equals(that.regionWithLowerCoordOnContig)) return false;
        if (!regionWithHigherCoordOnContig.equals(that.regionWithHigherCoordOnContig)) return false;
        if (strandSwitch != that.strandSwitch) return false;
        if (!insertionMappings.equals(that.insertionMappings)) return false;
        return goodNonCanonicalMappingSATag.equals(that.goodNonCanonicalMappingSATag);
    }

    @Override
    public int hashCode() {
        int result = sourceContigName.hashCode();
        result = 31 * result + regionWithLowerCoordOnContig.hashCode();
        result = 31 * result + regionWithHigherCoordOnContig.hashCode();
        result = 31 * result + strandSwitch.ordinal();
        result = 31 * result + (isForwardStrandRepresentation ? 1 : 0);
        result = 31 * result + insertionMappings.hashCode();
        result = 31 * result + goodNonCanonicalMappingSATag.hashCode();
        return result;
    }
}
