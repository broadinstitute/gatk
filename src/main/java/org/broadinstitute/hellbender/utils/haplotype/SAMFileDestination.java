package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;

/**
 * Class used to direct output from a HaplotypeBAMWriter to a bam/sam file.
 */
public final class SAMFileDestination extends HaplotypeBAMDestination {
    private final SAMFileWriter samWriter;

    /**
     * Create a new file destination.
     *
     * @param sourceHeader SAMFileHeader used to seed the output SAMFileHeader for this destination, must not be null
     * @param haplotypeReadGroupID read group ID used when writing haplotypes as reads
     */
    public SAMFileDestination(File outFile, SAMFileHeader sourceHeader, String haplotypeReadGroupID) {
        super(sourceHeader, haplotypeReadGroupID);
        samWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(getBAMOutputHeader(), false, outFile);
    }

    /**
     * Close any resources associated with this destination.
     */
    @Override
    void close() { samWriter.close(); }

    /**
     * Write a read to the output file specified by this destination.
     *
     * @param read the read to write out, must not be null
     */
    @Override
    public void add(final GATKRead read) {
        Utils.nonNull(read, "read cannot be null");
        samWriter.addAlignment(read.convertToSAMRecord(getBAMOutputHeader()));
    }
}