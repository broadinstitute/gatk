package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Arrays;
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
class ChimericAlignment {

    static final String NO_SEQUENCE = "none";
    String contigId;
    AlignmentRegion region1;
    AlignmentRegion region2;
    String homology;
    String insertedSequence;
    List<String> insertionMappings;

    /**
     * Construct a new ChimericAlignment from two AlignmentRegions.
     * Assumes {@code region1} has a lower {@link AlignmentRegion#startInAssembledContig} than {@code region2}.
     */
    protected ChimericAlignment(final AlignmentRegion region1, final AlignmentRegion region2,
                                final String insertedSequence, final String homology,
                                final List<String> insertionMappings) {

        final String contigId = region1.contigId;
        Utils.validateArg(contigId.equals(region2.contigId), "two alignment regions used to construct chimeric alignment are not from the same assembled contig.");

        this.contigId = contigId;
        this.region1 = region1;
        this.region2 = region2;
        this.homology = homology;
        this.insertedSequence = insertedSequence;
        this.insertionMappings = insertionMappings;
    }

    @SuppressWarnings("unchecked")
    protected ChimericAlignment(final Kryo kryo, final Input input) {
        this.contigId = input.readString();
        this.region1 = kryo.readObject(input, AlignmentRegion.class);
        this.region2 = kryo.readObject(input, AlignmentRegion.class);
        this.homology = input.readString();
        this.insertedSequence = input.readString();
        this.insertionMappings = (ArrayList<String>) kryo.readObject(input, ArrayList.class);
    }

    /**
     * TODO: never used/tested method, should remove?
     *  Parses a tab-delimited assembled breakpoint line into an ChimericAlignment object. Fields should be in the same order as that produced by toString():
     *
     *  contigId
     *  alignmentRegion1.toString()
     *  alignmentRegion2.toString()
     *  insertedSequence
     *  homology
     *
     */
    protected ChimericAlignment fromString(String assembledBreakpointLine) {
        final String[] fields = assembledBreakpointLine.split("\t");
        return fromFields(fields);
    }

    // TODO: private method that's only called in fromString, which is never accessed. Should remove?
    private ChimericAlignment fromFields(final String[] fields) {
        try {
//            final String contigId = fields[0].replaceFirst("^>","");
            final String[] alignmentRegion1Fields = Arrays.copyOfRange(fields, 1, 10);
            final AlignmentRegion alignmentRegion1 = AlignmentRegion.fromString(alignmentRegion1Fields);
            final String[] alignmentRegion2Fields = Arrays.copyOfRange(fields, 10, 19);
            final AlignmentRegion alignmentRegion2 = AlignmentRegion.fromString(alignmentRegion2Fields);
            final String insertedSequence = fields[19].equals(NO_SEQUENCE) ? "" : fields[19];
            final String homology = fields[20].equals(NO_SEQUENCE) ? "" : fields[20];
            final List<String> insertionMappings = Arrays.asList(fields[21].split(";"));
            return new ChimericAlignment(alignmentRegion1, alignmentRegion2, homology, insertedSequence, insertionMappings);
        } catch (final NumberFormatException nfe) {
            throw new GATKException(Arrays.toString(fields), nfe);
        }
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(contigId);
        kryo.writeObject(output, region1);
        kryo.writeObject(output, region2);
        output.writeString(homology);
        output.writeString(insertedSequence);
        kryo.writeObject(output, insertionMappings);
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
    public String toString() {
        return contigId +
                "\t" +
                region1.toString() +
                "\t" +
                region2.toString() +
                "\t" +
                ("".equals(insertedSequence) ? NO_SEQUENCE : insertedSequence) +
                "\t" +
                ("".equals(homology) ? NO_SEQUENCE : homology);
    }

    /**
     * A naive way to test if a chimeric alignment supports an inversion event.
     * Returning true does not necessarily mean there's actually an inversion.
     */
    @VisibleForTesting
    boolean involvesStrandSwitch() {
        return region1.forwardStrand != region2.forwardStrand;
    }

    /**
     * Returns the reference coordinates of the left and right breakpoints implied by this chimeric alignment.
     * If there is homologous sequence represented in the alignments, it will be assigned to the side of the breakpoint
     * with higher reference coordinates.
     */
    @VisibleForTesting
    final Tuple2<SimpleInterval, SimpleInterval> getLeftAndRightBreakpointsOnReferenceLeftAlignedForHomology() {
        final int leftBreakpointCoord, rightBreakpointCoord;
        if (region1.referenceInterval.getStart() < region2.referenceInterval.getStart()) {
            leftBreakpointCoord = region1.forwardStrand ? region1.referenceInterval.getEnd() - homology.length() : region1.referenceInterval.getStart();
            rightBreakpointCoord = region2.forwardStrand ? region2.referenceInterval.getStart() + homology.length() : region2.referenceInterval.getEnd();
        } else {
            leftBreakpointCoord = region2.forwardStrand ? region2.referenceInterval.getStart() : region2.referenceInterval.getEnd() - homology.length();
            rightBreakpointCoord = region1.forwardStrand ? region1.referenceInterval.getEnd() : region1.referenceInterval.getStart() + homology.length();
        }

        final SimpleInterval leftBreakpoint = new SimpleInterval(region1.referenceInterval.getContig(), leftBreakpointCoord, leftBreakpointCoord);
        final SimpleInterval rightBreakpoint = new SimpleInterval(region2.referenceInterval.getContig(), rightBreakpointCoord, rightBreakpointCoord);
        return new Tuple2<>(leftBreakpoint, rightBreakpoint);
    }
}
