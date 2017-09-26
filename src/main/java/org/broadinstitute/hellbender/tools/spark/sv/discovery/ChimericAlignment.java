package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;

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
                             final List<String> insertionMappings, final String sourceContigName) {

        this.sourceContigName = sourceContigName;

        this.regionWithLowerCoordOnContig = intervalWithLowerCoordOnContig;
        this.regionWithHigherCoordOnContig = intervalWithHigherCoordOnContig;

        this.strandSwitch = determineStrandSwitch(intervalWithLowerCoordOnContig, intervalWithHigherCoordOnContig);

        final boolean involvesRefIntervalSwitch = involvesRefPositionSwitch(intervalWithLowerCoordOnContig, intervalWithHigherCoordOnContig);
        this.isForwardStrandRepresentation = isForwardStrandRepresentation(intervalWithLowerCoordOnContig, intervalWithHigherCoordOnContig,
                this.strandSwitch, involvesRefIntervalSwitch);

        this.insertionMappings = insertionMappings;
    }

    // TODO: 12/14/16 Skipping simple translocation events
    /**
     * Parse all alignment records for a single locally-assembled contig and generate chimeric alignments if available.
     * Applies certain filters to skip the input alignment regions that are:
     *     1) if the alignment region's mapping quality is below a certain threshold, it is skipped
     *     2) if the alignment region is too small, it is skipped
     * If the alignment region passes the above two filters and the next alignment region could be treated as potential inserted sequence,
     * note down the mapping & alignment information of that region and skip it
     * @param alignedContig                     made of (sorted {alignmentIntervals}, sequence) of a potentially-signalling locally-assembled contig
     * @param uniqueRefSpanThreshold            for an alignment interval to be used to construct a ChimericAlignment,
     * @param filterAlignmentByMqOrLength       whether to perform filtering based on MQ for first several alignments and based alignment length for all alignments
     * @param filterWhollyContainedAlignments   whether to perform filtering on alignments whose ref span is contained in others'
     */
    @VisibleForTesting
    public static List<ChimericAlignment> parseOneContig(final AlignedContig alignedContig,
                                                         final int uniqueRefSpanThreshold,
                                                         final boolean filterAlignmentByMqOrLength,
                                                         final boolean filterWhollyContainedAlignments) {

        if (alignedContig.alignmentIntervals.size() < 2) {
            return new ArrayList<>();
        }

        final Iterator<AlignmentInterval> iterator = alignedContig.alignmentIntervals.iterator();

        // fast forward to the first alignment region with high MapQ
        AlignmentInterval current = iterator.next();
        if (filterAlignmentByMqOrLength) {
            while (mapQualTooLow(current, CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD) && iterator.hasNext()) {
                current = iterator.next();
            }
        }

        final List<ChimericAlignment> results = new ArrayList<>(alignedContig.alignmentIntervals.size() - 1);
        final List<String> insertionMappings = new ArrayList<>();

        // then iterate over the AR's in pair to identify CA's.
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            if (filterAlignmentByMqOrLength && firstAlignmentIsTooShort(current, next, uniqueRefSpanThreshold)) {
                continue;
            } else if (nextAlignmentMayBeInsertion(current, next, uniqueRefSpanThreshold, CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, filterWhollyContainedAlignments)) {
                if (iterator.hasNext()) {
                    insertionMappings.add(next.toPackedString());
                    continue;
                } else {
                    break;
                }
            }

            final boolean isNotSimpleTranslocation = isNotSimpleTranslocation(current, next,
                    determineStrandSwitch(current, next), involvesRefPositionSwitch(current, next));
            if (isNotSimpleTranslocation)
                results.add(new ChimericAlignment(current, next, insertionMappings, alignedContig.contigName));

            current = next;
        }

        return results;
    }

    // TODO: 11/22/16 it might also be suitable to consider the reference context this alignment region is mapped to
    //       and not simply apply a hard filter (need to think about how to test)
    private static boolean mapQualTooLow(final AlignmentInterval next, final int mapQThresholdInclusive) {
        return next.mapQual < mapQThresholdInclusive;
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
                                                      final Integer minAlignLength, final int mapQThresholdInclusive,
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
     * Determine if the region that maps to a lower coordinate on the contig also maps to a lower coordinate on the reference.
     */
    @VisibleForTesting
    public static boolean involvesRefPositionSwitch(final AlignmentInterval regionWithLowerCoordOnContig,
                                                    final AlignmentInterval regionWithHigherCoordOnContig) {

        return regionWithHigherCoordOnContig.referenceSpan.getStart() < regionWithLowerCoordOnContig.referenceSpan.getStart();
    }

    // TODO: 12/15/16 simple translocations are defined here and at this time as inter-chromosomal ones and
    //                intra-chromosomal translocations that involve strand switch
    //                to get the 2nd case right, we need evidence flanking both sides of the inversion, and that could be difficult
    // TODO: 1/17/17 this is used for filtering out possible translocations that we are not currently handling,
    //                but it overkills some insertions where the inserted sequence maps to chromosomes other than that of the flanking regions
    /**
     * Determine if the two regions indicate a inter-chromosomal translocation or intra-chromosomal translocation that
     * DOES NOT involve a strand switch.
     */
    @VisibleForTesting
    public static boolean isNotSimpleTranslocation(final AlignmentInterval regionWithLowerCoordOnContig,
                                                   final AlignmentInterval regionWithHigherCoordOnContig,
                                                   final StrandSwitch strandSwitch,
                                                   final boolean involvesReferenceIntervalSwitch) {
        // logic is: must be the same reference chromosome for it not to be an inter-chromosomal translocation
        //           and when regions are mapped to the same reference chromosome, there cannot be reference position swap
        final boolean sameChromosome = regionWithLowerCoordOnContig.referenceSpan.getContig()
                                       .equals(regionWithHigherCoordOnContig.referenceSpan.getContig());
        return sameChromosome
                &&
                (strandSwitch!=StrandSwitch.NO_SWITCH
                        || involvesReferenceIntervalSwitch == !regionWithLowerCoordOnContig.forwardStrand);
    }

    /**
     * An SV event could be detected from a contig that seem originate from the forward or reverse strand of the reference,
     * besides the annotation that the alignment flanking regions might flank either side of the two breakpoints.
     */
    @VisibleForTesting
    static boolean isForwardStrandRepresentation(final AlignmentInterval regionWithLowerCoordOnContig,
                                                 final AlignmentInterval regionWithHigherCoordOnContig,
                                                 final StrandSwitch strandSwitch,
                                                 final boolean involvesReferenceIntervalSwitch) {

        if (strandSwitch == StrandSwitch.NO_SWITCH) {
            return regionWithLowerCoordOnContig.forwardStrand && regionWithHigherCoordOnContig.forwardStrand;
        } else {
            return !involvesReferenceIntervalSwitch;
        }
    }

    public boolean isNotSimpleTranslocation() {
        return isNotSimpleTranslocation(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig, strandSwitch,
                                        involvesRefPositionSwitch(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig));
    }

    Tuple2<SimpleInterval, SimpleInterval> getCoordSortedReferenceSpans() {
        if (involvesRefPositionSwitch(regionWithLowerCoordOnContig, regionWithHigherCoordOnContig))
            return new Tuple2<>(regionWithHigherCoordOnContig.referenceSpan, regionWithLowerCoordOnContig.referenceSpan);
        else
            return new Tuple2<>(regionWithLowerCoordOnContig.referenceSpan, regionWithHigherCoordOnContig.referenceSpan);
    }


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
