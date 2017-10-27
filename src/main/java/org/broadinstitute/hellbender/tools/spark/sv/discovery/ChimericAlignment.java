package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

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

    final String sourceContigName;

    final AlignmentInterval regionWithLowerCoordOnContig;
    final AlignmentInterval regionWithHigherCoordOnContig;

    final StrandSwitch strandSwitch;
    final boolean isForwardStrandRepresentation;

    final List<String> insertionMappings;

    public List<AlignmentInterval> getAlignmentIntervals() {
        return Arrays.asList(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig);
    }

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
    //////////// BELOW ARE CODE PATH USED FOR INSERTION, DELETION, AND DUPLICATION (INV OR NOT) AND INVERSION, AND ARE TESTED ONLY FOR THAT PURPOSE

    /**
     * Parse all alignment records for a single locally-assembled contig and generate chimeric alignments if available.
     * Applies certain filters to skip the input alignment regions that are:
     *     1) if the alignment region's mapping quality is below a certain threshold, it is skipped
     *     2) if the alignment region is too small, it is skipped
     * If the alignment region passes the above two filters and the next alignment region could be treated as potential inserted sequence,
     * note down the mapping & alignment information of that region and skip it
     * @param alignedContig          made of (sorted {alignmentIntervals}, sequence) of a potentially-signalling locally-assembled contig
     * @param referenceDictionary    reference sequence dictionary
     * @param filterAlignmentByMqOrLength
     * @param uniqueRefSpanThreshold for an alignment interval to be used to construct a ChimericAlignment,
*                               how long a unique--i.e. the same ref span is not covered by other alignment intervals--alignment on the reference must it have
     * @param mapQualThresholdInclusive
     * @param filterWhollyContainedAlignments
     */
    @VisibleForTesting
    public static List<ChimericAlignment> parseOneContig(final AlignedContig alignedContig,
                                                         final SAMSequenceDictionary referenceDictionary,
                                                         final boolean filterAlignmentByMqOrLength,
                                                         final int uniqueRefSpanThreshold,
                                                         final int mapQualThresholdInclusive,
                                                         final boolean filterWhollyContainedAlignments) {

        if (alignedContig.alignmentIntervals.size() < 2) {
            return new ArrayList<>();
        }

        final Iterator<AlignmentInterval> iterator = alignedContig.alignmentIntervals.iterator();

        // fast forward to the first alignment region with high MapQ
        AlignmentInterval current = iterator.next();
        if (filterAlignmentByMqOrLength) {
            while (mapQualTooLow(current, mapQualThresholdInclusive) && iterator.hasNext()) {
                current = iterator.next();
            }
        }

        final List<ChimericAlignment> results = new ArrayList<>(alignedContig.alignmentIntervals.size() - 1);
        final List<String> insertionMappings = new ArrayList<>();

        // then iterate over the AR's in pair to identify CA's.
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            if (filterAlignmentByMqOrLength) {
                if (firstAlignmentIsTooShort(current, next, uniqueRefSpanThreshold)) {
                    continue;
                } else if (nextAlignmentMayBeInsertion(current, next, mapQualThresholdInclusive, uniqueRefSpanThreshold, filterWhollyContainedAlignments)) {
                    if (iterator.hasNext()) {
                        insertionMappings.add(next.toPackedString());
                        continue;
                    } else {
                        break;
                    }
                }
            }

            // TODO: 10/18/17 this way of filtering CA based on not quality but alignment characteristics is temporary:
            //       this was initially developed for ins/del (and tested for that purpose), simple translocations travel through a different code path at the moment.
            // TODO: ultimately we need to merge these two code paths
            final ChimericAlignment chimericAlignment = new ChimericAlignment(current, next, insertionMappings,
                    alignedContig.contigName, referenceDictionary);
            // the following check/filter is due to the fact that simple translocations are to be handled in a different code path
            if (chimericAlignment.isNotSimpleTranslocation())
                results.add(chimericAlignment);

            current = next;
        }

        return results;
    }

    //==================================================================================================================
    // TODO: 11/22/16 it might also be suitable to consider the reference context this alignment region is mapped to
    //       and not simply apply a hard filter (need to think about how to test)
    private static boolean mapQualTooLow(final AlignmentInterval aln, final int mapQThresholdInclusive) {
        return aln.mapQual < mapQThresholdInclusive;
    }

    /**
     * @return if {@code first} is too short, when considering overlap with {@code second}
     */
    @VisibleForTesting
    static boolean firstAlignmentIsTooShort(final AlignmentInterval first, final AlignmentInterval second,
                                            final Integer minAlignLength) {
        return first.referenceSpan.size() - AlignmentInterval.overlapOnContig(first, second) < minAlignLength;
    }

    /**
     * To implement the idea that for two consecutive alignment regions of a contig, the one with higher reference coordinate might be a novel insertion.
     */
    @VisibleForTesting
    public static boolean nextAlignmentMayBeInsertion(final AlignmentInterval current, final AlignmentInterval next,
                                                      final Integer mapQThresholdInclusive, final Integer minAlignLength,
                                                      final boolean filterWhollyContained) {
        // not unique: inserted sequence may have low mapping quality (low reference uniqueness) or may be very small (low read uniqueness)
        final boolean isNotUnique = mapQualTooLow(next, mapQThresholdInclusive) || firstAlignmentIsTooShort(next, current, minAlignLength);
        return isNotUnique
                ||
                (filterWhollyContained && (current.referenceSpan.contains(next.referenceSpan) || next.referenceSpan.contains(current.referenceSpan)));
    }

    @VisibleForTesting
    public static StrandSwitch determineStrandSwitch(final AlignmentInterval first, final AlignmentInterval second) {
        if (first.forwardStrand == second.forwardStrand) {
            return StrandSwitch.NO_SWITCH;
        } else {
            return first.forwardStrand ? StrandSwitch.FORWARD_TO_REVERSE : StrandSwitch.REVERSE_TO_FORWARD;
        }
    }

    /**
     * An SV event could be detected from a contig that seem originate from the forward or reverse strand of the reference,
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

    /**
     * See {@link #isNotSimpleTranslocation()} for logic.
     */
    @VisibleForTesting
    static boolean isLikelySimpleTranslocation(final AlignmentInterval regionWithLowerCoordOnContig,
                                               final AlignmentInterval regionWithHigherCoordOnContig,
                                               final StrandSwitch strandSwitch) {

        if (!regionWithLowerCoordOnContig.referenceSpan.getContig()
                .equals(regionWithHigherCoordOnContig.referenceSpan.getContig()))
            return true;

        if ( !strandSwitch.equals(StrandSwitch.NO_SWITCH) )
            return false;

        final SimpleInterval referenceSpanOne = regionWithLowerCoordOnContig.referenceSpan,
                             referenceSpanTwo = regionWithHigherCoordOnContig.referenceSpan;

        if (referenceSpanOne.contains(referenceSpanTwo) || referenceSpanTwo.contains(referenceSpanOne))
            return false;

        if (regionWithLowerCoordOnContig.forwardStrand) {
            return referenceSpanOne.getStart() > referenceSpanTwo.getEnd();
        } else {
            return referenceSpanTwo.getStart() > referenceSpanOne.getEnd();
        }
    }

    /**
     * Determine if the chimeric alignment indicates a simple translocation.
     * Simple translocations are defined here and at this time as:
     * <ul>
     *  <li>inter-chromosomal translocations, i.e. novel adjacency between different reference chromosomes, or</li>
     *  <li>intra-chromosomal translocation that DOES NOT involve a strand switch, i.e.
     *      novel adjacency between reference locations on the same chromosome involving NO strand switch,
     *      but in the meantime, the two inducing alignments CANNOT overlap each other since that would point to
     *      incomplete picture, hence not "simple" anymore.
     *  </li>
     * </ul>
     * A caveat is that this does not cover the case when the novel adjacency suggested by the CA is between
     * two reference locations on the same chromosome, but involves a strand switch,
     * which could be a translocation or inversion breakpoint.
     * But to fully resolve this case, we need other types of evidence, hence should not be the task of this function.
     */
    public boolean isNotSimpleTranslocation() {
        return !isLikelySimpleTranslocation(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig, strandSwitch);
    }

    /**
     * See {@link #isLikelyInvertedDuplication()}
     */
    @VisibleForTesting
    public static boolean isLikelyInvertedDuplication(final AlignmentInterval one, final AlignmentInterval two) {
        if (one.forwardStrand == two.forwardStrand)
            return false;
        return 2 * AlignmentInterval.overlapOnRefSpan(one, two) >
                Math.min(one.endInAssembledContig - one.startInAssembledContig,
                        two.endInAssembledContig - two.startInAssembledContig) + 1;
    }

    /**
     * todo : see ticket #3529
     * @return true iff the two AI of the {@code longRead} are
     *         1) of different strand and
     *         2) overlap on reference is more than half of the two AI's minimal read span.
     */
    public boolean isLikelyInvertedDuplication() {
        return isLikelyInvertedDuplication(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig);
    }

    //==================================================================================================================
    public String onErrStringRep() {
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
        result = 31 * result + strandSwitch.hashCode();
        result = 31 * result + (isForwardStrandRepresentation ? 1 : 0);
        result = 31 * result + insertionMappings.hashCode();
        return result;
    }
}
