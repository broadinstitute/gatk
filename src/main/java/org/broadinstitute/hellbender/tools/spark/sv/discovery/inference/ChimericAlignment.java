package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Holds information about a split alignment of a contig, which may represent an SV breakpoint. Each ChimericAlignment
 * represents the junction on the contig of two aligned regions. For example, if a contig aligns to three different regions
 * of the genome (with one primary and two supplementary alignment records), there will be two ChimericAlignment
 * objects created, one to represent each junction between alignment regions:
 *
 * Example Contig:
 * ACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTG
 * Alignment regions:
 * |---------1:100-200------------|
 *                                 |----------2:100-200------------------|
 *                                                                       |----------3:100-200-----------------|
 * Assembled breakpoints:
 * 1) links 1:100-200 to 2:100-200
 * 2) links 2:100-200 to 3:100-200
 *
 * Inserted sequence contains portions of the contig that are aligned to neither region, and therefore may be inserted in
 * the sample. For example, a translocation breakpoint with a micro-insertion:
 *
 * Contig:
 * ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
 * Alignment regions:
 * |-----1:100-200-------|
 *                          |----2:100-200-----|
 * Inserted sequence:
 *  GA
 *
 * Homology represents ambiguity about the exact location of the breakpoint. For example, in this case one alignment
 * region ends with "AC" and the next begins with AC, so we don't know if the AC truly belongs with the first or
 * second alignment region.
 *
 * Contig:
 * ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
 * Alignment regions:
 * |-----1:100-200-------|
 *                    |-----2:100-200----------|
 * Homology:
 *  AC
 */
@DefaultSerializer(ChimericAlignment.Serializer.class)
public class ChimericAlignment {

    public final String sourceContigName;

    public final AlignmentInterval regionWithLowerCoordOnContig;
    public final AlignmentInterval regionWithHigherCoordOnContig;

    final StrandSwitch strandSwitch;
    final boolean isForwardStrandRepresentation;

    public final List<String> insertionMappings;

    protected ChimericAlignment(final Kryo kryo, final Input input) {

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
    }

    /**
     * Construct a new ChimericAlignment from two alignment intervals.
     * Assumes {@code intervalWithLowerCoordOnContig} has a lower {@link AlignmentInterval#startInAssembledContig}
     * than {@code regionWithHigherCoordOnContig}.
     */
    @VisibleForTesting
    public ChimericAlignment(final AlignmentInterval intervalWithLowerCoordOnContig, final AlignmentInterval intervalWithHigherCoordOnContig,
                             final List<String> insertionMappings, final String sourceContigName,
                             final SAMSequenceDictionary referenceDictionary) {

        this.sourceContigName = sourceContigName;

        this.regionWithLowerCoordOnContig = intervalWithLowerCoordOnContig;
        this.regionWithHigherCoordOnContig = intervalWithHigherCoordOnContig;

        this.strandSwitch = determineStrandSwitch(intervalWithLowerCoordOnContig, intervalWithHigherCoordOnContig);

        this.isForwardStrandRepresentation =
                isForwardStrandRepresentation(intervalWithLowerCoordOnContig, intervalWithHigherCoordOnContig, strandSwitch, referenceDictionary);

        this.insertionMappings = insertionMappings;
    }

    // =================================================================================================================

    /**
     * See
     * {@link AssemblyContigAlignmentSignatureClassifier#isCandidateSimpleTranslocation(AlignmentInterval, AlignmentInterval, StrandSwitch)}
     * for definition of "simple translocations".
     */
    boolean isLikelySimpleTranslocation() {
        return AssemblyContigAlignmentSignatureClassifier.isCandidateSimpleTranslocation(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig, strandSwitch);
    }

    /**
     * See {@link AssemblyContigAlignmentSignatureClassifier#isCandidateInvertedDuplication(AlignmentInterval, AlignmentInterval)}
     */
    boolean isLikelyInvertedDuplication() {
        return AssemblyContigAlignmentSignatureClassifier.isCandidateInvertedDuplication(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig);
    }

    public boolean isNeitherSimpleTranslocationNorIncompletePicture() {
        return !isLikelySimpleTranslocation()
                &&
                !AssemblyContigAlignmentSignatureClassifier.hasIncompletePictureFromTwoAlignments(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig);
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
     * {@link #isNeitherSimpleTranslocationNorIncompletePicture()} and
     * {@link #determineStrandSwitch(AlignmentInterval, AlignmentInterval)} == {@link StrandSwitch#NO_SWITCH}
     */
    static final class DistancesBetweenAlignmentsOnRefAndOnRead {
        final int distBetweenAlignRegionsOnRef; // distance-1 between the two regions on reference, denoted as d1 in the comments below
        final int distBetweenAlignRegionsOnCtg; // distance-1 between the two regions on contig, denoted as d2 in the comments below

        final int leftAlnRefEnd;
        final int rightAlnRefStart;
        final int firstAlnCtgEnd;
        final int secondAlnCtgStart;

        DistancesBetweenAlignmentsOnRefAndOnRead(final int distBetweenAlignRegionsOnRef,
                                                 final int distBetweenAlignRegionsOnCtg,
                                                 final int leftAlnRefEnd,
                                                 final int rightAlnRefStart,
                                                 final int firstAlnCtgEnd,
                                                 final int secondAlnCtgStart) {
            this.distBetweenAlignRegionsOnRef = distBetweenAlignRegionsOnRef;
            this.distBetweenAlignRegionsOnCtg = distBetweenAlignRegionsOnCtg;
            this.leftAlnRefEnd = leftAlnRefEnd;
            this.rightAlnRefStart = rightAlnRefStart;
            this.firstAlnCtgEnd = firstAlnCtgEnd;
            this.secondAlnCtgStart = secondAlnCtgStart;
        }
    }

    DistancesBetweenAlignmentsOnRefAndOnRead getDistancesBetweenAlignmentsOnRefAndOnRead() {
        Utils.validateArg(isNeitherSimpleTranslocationNorIncompletePicture(),
                "Assumption that the simple chimera is neither incomplete picture nor simple translocation is violated.\n" +
                        toString());
        Utils.validateArg(strandSwitch.equals(StrandSwitch.NO_SWITCH),
                "Assumption that the simple chimera is neither incomplete picture nor simple translocation is violated.\n" +
                        toString());
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

    /**
     * Roughly similar to {@link DiscoverVariantsFromContigAlignmentsSAMSpark#nextAlignmentMayBeInsertion(AlignmentInterval, AlignmentInterval, Integer, Integer, boolean)}:
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

    @Override
    public String toString() {
        return sourceContigName +
                "\t" +
                regionWithLowerCoordOnContig.toPackedString() +
                "\t" +
                regionWithHigherCoordOnContig.toPackedString();
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
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ChimericAlignment> {
        @Override
        public void write(final Kryo kryo, final Output output, final ChimericAlignment chimericAlignment) {
            chimericAlignment.serialize(kryo, output);
        }

        @Override
        public ChimericAlignment read(final Kryo kryo, final Input input, final Class<ChimericAlignment> klass ) {
            return new ChimericAlignment(kryo, input);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ChimericAlignment that = (ChimericAlignment) o;

        if (isForwardStrandRepresentation != that.isForwardStrandRepresentation) return false;
        if (!sourceContigName.equals(that.sourceContigName)) return false;
        if (!regionWithLowerCoordOnContig.equals(that.regionWithLowerCoordOnContig)) return false;
        if (!regionWithHigherCoordOnContig.equals(that.regionWithHigherCoordOnContig)) return false;
        if (strandSwitch != that.strandSwitch) return false;
        return insertionMappings.equals(that.insertionMappings);
    }

    @Override
    public int hashCode() {
        int result = sourceContigName.hashCode();
        result = 31 * result + regionWithLowerCoordOnContig.hashCode();
        result = 31 * result + regionWithHigherCoordOnContig.hashCode();
        result = 31 * result + strandSwitch.ordinal();
        result = 31 * result + (isForwardStrandRepresentation ? 1 : 0);
        result = 31 * result + insertionMappings.hashCode();
        return result;
    }
}
