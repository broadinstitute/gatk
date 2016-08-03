package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Arrays;
import java.util.List;

/**
 * Holds information about a split alignment of a contig, which may represent an SV breakpoint. Each AssembledBreakpoint
 * represents the junction on the contig of two aligned regions. For example, if a contig aligns to three different regions
 * of the genome (with one primary and two supplementary alignment records), there will be two AssembledBreakpoint
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
class AssembledBreakpoint {
    String contigId;
    ContigAligner.AlignmentRegion region1;
    ContigAligner.AlignmentRegion region2;
    String insertedSequence;
    String homology;
    List<String> insertionMappings;

    public AssembledBreakpoint(final String contigId, final ContigAligner.AlignmentRegion region1, final ContigAligner.AlignmentRegion region2, final String insertedSequence, final String homology, final List<String> insertionMappings) {
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
                ("".equals(insertedSequence) ? "NA" : insertedSequence) +
                "\t" +
                ("".equals(homology) ? "NA" : homology);
    }

    /**
     *  Parses a tab-delimited assembled breakpoint line into an AssembledBreakpoint object
     */
    public static AssembledBreakpoint fromString(String assembledBreakpointLine) {
        final String[] fields = assembledBreakpointLine.split("\t");
        return fromFields(fields);
    }

    public static AssembledBreakpoint fromFields(final String[] fields) {
        try {
            final String contigId = fields[0].replaceFirst("^>","");
            final String[] alignmentRegion1Fields = Arrays.copyOfRange(fields, 1, 10);
            final ContigAligner.AlignmentRegion alignmentRegion1 = ContigAligner.AlignmentRegion.fromString(alignmentRegion1Fields);
            final String[] alignmentRegion2Fields = Arrays.copyOfRange(fields, 10, 19);
            final ContigAligner.AlignmentRegion alignmentRegion2 = ContigAligner.AlignmentRegion.fromString(alignmentRegion2Fields);
            final String insertedSequence = fields[19].equals("NA") ? "" : fields[19];
            final String homology = fields[20].equals("NA") ? "" : fields[20];
            final List<String> insertionMappings = Arrays.asList(fields[21].split(";"));
            return new AssembledBreakpoint(contigId, alignmentRegion1, alignmentRegion2, insertedSequence, homology, insertionMappings);
        } catch (final NumberFormatException nfe) {
            throw new GATKException(Arrays.toString(fields), nfe);
        }

    }

    public SimpleInterval getLeftAlignedLeftBreakpointOnAssembledContig() {
        final int alignmentEnd = region1.forwardStrand ? region1.referenceInterval.getEnd() : region1.referenceInterval.getStart();
        final int position = region1.forwardStrand ? alignmentEnd - homology.length() : alignmentEnd;
        return new SimpleInterval(region1.referenceInterval.getContig(), position, position);
    }

    public SimpleInterval getLeftAlignedRightBreakpointOnAssembledContig() {
        final int alignmentStart = region2.forwardStrand ? region2.referenceInterval.getStart() : region2.referenceInterval.getEnd();
        final int position = region2.forwardStrand ? alignmentStart : alignmentStart - homology.length();
        return new SimpleInterval(region2.referenceInterval.getContig(), position, position);
    }

    public ContigAligner.BreakpointAllele getBreakpointAllele() {
        final SimpleInterval leftAlignedLeftBreakpointOnAssembledContig = getLeftAlignedLeftBreakpointOnAssembledContig();
        final SimpleInterval leftAlignedRightBreakpointOnAssembledContig = getLeftAlignedRightBreakpointOnAssembledContig();

        final boolean isFiveToThreeInversion = region1.forwardStrand && ! region2.forwardStrand &&
                ! region1.referenceInterval.contains(region2.referenceInterval);

        final boolean isThreeToFiveInversion = ! region1.forwardStrand && region2.forwardStrand
                && ! region1.referenceInterval.contains(region2.referenceInterval);

        if (! leftAlignedLeftBreakpointOnAssembledContig.getContig().equals(leftAlignedRightBreakpointOnAssembledContig.getContig())) {
            return new ContigAligner.BreakpointAllele(leftAlignedLeftBreakpointOnAssembledContig, leftAlignedRightBreakpointOnAssembledContig, insertedSequence, homology, isFiveToThreeInversion, isThreeToFiveInversion);
        } else if ( leftAlignedLeftBreakpointOnAssembledContig.getStart() < leftAlignedRightBreakpointOnAssembledContig.getStart()) {
            return new ContigAligner.BreakpointAllele(leftAlignedLeftBreakpointOnAssembledContig, leftAlignedRightBreakpointOnAssembledContig, insertedSequence, homology, isFiveToThreeInversion, isThreeToFiveInversion);
        } else {
            return new ContigAligner.BreakpointAllele(leftAlignedRightBreakpointOnAssembledContig, leftAlignedLeftBreakpointOnAssembledContig, insertedSequence, homology, isFiveToThreeInversion, isThreeToFiveInversion);
        }
    }
}
