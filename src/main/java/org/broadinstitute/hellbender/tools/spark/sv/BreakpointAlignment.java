package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.util.Arrays;
import java.util.List;

/**
 * Holds information about a split alignment of a contig, which may represent an SV breakpoint. Each BreakpointAlignment
 * represents the junction on the contig of two aligned regions. For example, if a contig aligns to three different regions
 * of the genome (with one primary and two supplementary alignment records), there will be two BreakpointAlignment
 * objects created, one to represent each junction between alignment regions:
 *
 * Example Contig:
 * ACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTGCACTGACTG
 * Alignment regions:
 * |---------1:100-200------------|
 *                                 |----------2:100-200------------------|
 *                                                                       |----------3:100-200-----------------|
 * Assmbled breakpoints:
 * 1) links 1:100-200 to 2:100-200
 * 2) links 2:100-200 to 3:100-200
 *
 * Inserted sequence contains portions of the contig that are aligned to neither region, and therefore may be inserted in
 * the sample. For example, a translocation breakpoint with a microinsertion:
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
class BreakpointAlignment {
    public static final String NO_SEQUENCE = "none";
    String contigId;
    AlignmentRegion region1;
    AlignmentRegion region2;
    String insertedSequence;
    String homology;
    List<String> insertionMappings;

    /**
     * Construct a new BreakpointAlignment from two AlignmentRegions.
     * Assumes {@code region1} has a lower {@link AlignmentRegion#startInAssembledContig} than {@code region2}.
     */
    public BreakpointAlignment(final String contigId, final AlignmentRegion region1, final AlignmentRegion region2, final String insertedSequence, final String homology, final List<String> insertionMappings) {
        this.contigId = contigId;
        this.region1 = region1;
        this.region2 = region2;
        this.insertedSequence = insertedSequence;
        this.homology = homology;
        this.insertionMappings = insertionMappings;
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
     *  Parses a tab-delimited assembled breakpoint line into an BreakpointAlignment object. Fields should be in the same order as that produced by toString():
     *
     *  contigId
     *  alignmentRegion1.toString()
     *  alignmentRegion2.toString()
     *  insertedSequence
     *  homology
     *
     */
    public static BreakpointAlignment fromString(String assembledBreakpointLine) {
        final String[] fields = assembledBreakpointLine.split("\t");
        return fromFields(fields);
    }

    private static BreakpointAlignment fromFields(final String[] fields) {
        try {
            final String contigId = fields[0].replaceFirst("^>","");
            final String[] alignmentRegion1Fields = Arrays.copyOfRange(fields, 1, 10);
            final AlignmentRegion alignmentRegion1 = AlignmentRegion.fromString(alignmentRegion1Fields);
            final String[] alignmentRegion2Fields = Arrays.copyOfRange(fields, 10, 19);
            final AlignmentRegion alignmentRegion2 = AlignmentRegion.fromString(alignmentRegion2Fields);
            final String insertedSequence = fields[19].equals(NO_SEQUENCE) ? "" : fields[19];
            final String homology = fields[20].equals(NO_SEQUENCE) ? "" : fields[20];
            final List<String> insertionMappings = Arrays.asList(fields[21].split(";"));
            return new BreakpointAlignment(contigId, alignmentRegion1, alignmentRegion2, insertedSequence, homology, insertionMappings);
        } catch (final NumberFormatException nfe) {
            throw new GATKException(Arrays.toString(fields), nfe);
        }

    }

    private AlignmentRegion getLeftAlignmentRegion() {
        if (region1.referenceInterval.getStart() < region2.referenceInterval.getStart()) {
            return region1;
        } else {
            return region2;
        }
    }

    /**
     * Returns the reference coordinates of the left and right breakpoints implied by this alignment. If there is homologous sequence
     * represented in the alignments, it will be assigned to the side of the breakpoint with higher reference coordinates.
     * @return
     */
    @VisibleForTesting
    Tuple2<SimpleInterval, SimpleInterval> getLeftAndRightBreakpointsOnReferenceLeftAlignedForHomology() {
        if (region1 == getLeftAlignmentRegion()) {
            final int leftBreakpointCoord = region1.forwardStrand ? region1.referenceInterval.getEnd() - homology.length() : region1.referenceInterval.getStart();
            final SimpleInterval leftBreakpoint = new SimpleInterval(region1.referenceInterval.getContig(), leftBreakpointCoord, leftBreakpointCoord);
            final int rightBreakpointCoord = region2.forwardStrand ? region2.referenceInterval.getStart() + homology.length() : region2.referenceInterval.getEnd();
            final SimpleInterval rightBreakpoint = new SimpleInterval(region2.referenceInterval.getContig(), rightBreakpointCoord, rightBreakpointCoord);
            return new Tuple2<>(leftBreakpoint, rightBreakpoint);
        } else {
            final int leftBreakpointCoord = region2.forwardStrand ? region2.referenceInterval.getStart() : region2.referenceInterval.getEnd() - homology.length();
            final SimpleInterval leftBreakpoint = new SimpleInterval(region1.referenceInterval.getContig(), leftBreakpointCoord, leftBreakpointCoord);
            final int rightBreakpointCoord = region1.forwardStrand ? region1.referenceInterval.getEnd() : region1.referenceInterval.getStart() + homology.length();
            final SimpleInterval rightBreakpoint = new SimpleInterval(region2.referenceInterval.getContig(), rightBreakpointCoord, rightBreakpointCoord);
            return new Tuple2<>(leftBreakpoint, rightBreakpoint);
        }
    }

    
    /**
     * Returns the canonical representation of the breakpoint implied by this split contig alignment,
     * including whether it is a 3-5 or 5-3 inversion, and the homology and inserted sequence at the
     * breakpoint. The two intervals returned are 1bp intervals indicating the exact breakpoint
     * location. If there is homology at the breakpoint, the breakpoint locations will be left
     * aligned.
     * @return
     */
    public BreakpointAllele getBreakpointAllele() {

        final Tuple2<SimpleInterval, SimpleInterval> leftAndRightBreakpointsOnReferenceLeftAlignedForHomology = getLeftAndRightBreakpointsOnReferenceLeftAlignedForHomology();
        final SimpleInterval leftAlignedLeftBreakpointOnAssembledContig = leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._1();
        final SimpleInterval leftAlignedRightBreakpointOnAssembledContig = leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._2();

        final boolean isFiveToThreeInversion = region1.forwardStrand && !region2.forwardStrand;
        final boolean isThreeToFiveInversion = !region1.forwardStrand && region2.forwardStrand;

        if (!leftAlignedLeftBreakpointOnAssembledContig.getContig().equals(leftAlignedRightBreakpointOnAssembledContig.getContig())) {
            return new BreakpointAllele(leftAlignedLeftBreakpointOnAssembledContig, leftAlignedRightBreakpointOnAssembledContig, insertedSequence, homology, isFiveToThreeInversion, isThreeToFiveInversion, insertionMappings);
        } else if (leftAlignedLeftBreakpointOnAssembledContig.getStart() < leftAlignedRightBreakpointOnAssembledContig.getStart()) {
            return new BreakpointAllele(leftAlignedLeftBreakpointOnAssembledContig, leftAlignedRightBreakpointOnAssembledContig, insertedSequence, homology, isFiveToThreeInversion, isThreeToFiveInversion, insertionMappings);
        } else {
            return new BreakpointAllele(leftAlignedRightBreakpointOnAssembledContig, leftAlignedLeftBreakpointOnAssembledContig, insertedSequence, homology, isFiveToThreeInversion, isThreeToFiveInversion, insertionMappings);
        }
    }
}
