package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import java.nio.file.Path;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;


import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Class used to direct output from a HaplotypeBAMWriter to a bam/sam file.
 *
 */

public class HaplotypeBAMDestination {
    private final SAMFileGATKReadWriter samWriter;
    private final SAMFileHeader bamOutputHeader;
    private final String haplotypeReadGroupID;
    private final static String haplotypeSampleTag = "HC";

    /**
     * Create a new HaplotypeBAMDestination
     *
     * @param outPath path where output is written (doesn't have to be local)
     * @param createBamOutIndex true to create an index file for the bamout
     * @param createBamOutMD5 true to create an md5 file for the bamout
     * @param sourceHeader SAMFileHeader used to seed the output SAMFileHeader for this destination, must not be null
     * @param haplotypeReadGroupID  read group ID used when writing haplotypes as reads
     */
    protected HaplotypeBAMDestination(
            final Path outPath,
            final boolean createBamOutIndex,
            final boolean createBamOutMD5,
            final SAMFileHeader sourceHeader,
            final String haplotypeReadGroupID)
    {
        Utils.nonNull(outPath, "outputPath cannot be null");
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

        samWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(
                outPath,
                null,
                getBAMOutputHeader(), // use the header derived from the source header by HaplotypeBAMDestination
                false,
                createBamOutIndex,
                createBamOutMD5
        ));
    }

    /**
     * Write a read to the output specified by this destination.
     *
     * @param read the read to write out
     */
    public void add(final GATKRead read){
        Utils.nonNull(read, "read cannot be null");
        samWriter.addRead(read);
    }

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
     * Close any resources associated with this destination.
     */

    void close(){
        samWriter.close();
    };

    /**
     * Get the SAMFileHeader that is used for writing the output for this destination.
     * @return output SAMFileHeader
     */
    public SAMFileHeader getBAMOutputHeader() {
        return bamOutputHeader;
    }

}
