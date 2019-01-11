package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

/**
 * Utility class that allows easy creation of destinations for the HaplotypeBAMWriters
 *
 */
public abstract class HaplotypeBAMDestination {
    private final SAMFileHeader bamOutputHeader;
    private final String haplotypeReadGroupID;
    private final static String haplotypeSampleTag = "HC";

    /**
     * Create a new HaplotypeBAMDestination
     *
     * @param sourceHeader SAMFileHeader used to seed the output SAMFileHeader for this destination.
     * @param haplotypeReadGroupID read group ID used when writing haplotypes as reads
     */
    protected HaplotypeBAMDestination(SAMFileHeader sourceHeader, final String haplotypeReadGroupID) {
        Utils.nonNull(sourceHeader, "sourceHeader cannot be null");
        Utils.nonNull(haplotypeReadGroupID, "haplotypeReadGroupID cannot be null");
        this.haplotypeReadGroupID = haplotypeReadGroupID;

        bamOutputHeader = new SAMFileHeader();
        bamOutputHeader.setSequenceDictionary(sourceHeader.getSequenceDictionary());
        bamOutputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        final List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        readGroups.addAll(sourceHeader.getReadGroups()); // include the original read groups

        // plus an artificial read group for the haplotypes
        final SAMReadGroupRecord rgRec = new SAMReadGroupRecord(getHaplotypeReadGroupID());
        rgRec.setSample(haplotypeSampleTag);
        rgRec.setSequencingCenter("BI");
        readGroups.add(rgRec);
        bamOutputHeader.setReadGroups(readGroups);

        bamOutputHeader.addProgramRecord(new SAMProgramRecord("HalpotypeBAMWriter"));
    }

    /**
     * Write a read to the output specified by this destination.
     *
     * @param read the read to write out
     */
    public abstract void add(final GATKRead read);

    /**
     * Get the read group ID that is used by this writer when writing halpotypes as reads.
     *
     * @return read group ID
     */
    public String getHaplotypeReadGroupID() { return haplotypeReadGroupID; }

    /**
     * Get the sample tag that is used by this writer when writing halpotypes as reads.
     *
     * @return sample tag
     */
    public String getHaplotypeSampleTag() { return haplotypeSampleTag; }

    /**
     * Close the destination
     *
     */
    abstract void close();

    /**
     * Get the SAMFileHeader that is used for writing the output for this destination.
     * @return output SAMFileHeader
     */
    public SAMFileHeader getBAMOutputHeader() {
        return bamOutputHeader;
    }

}