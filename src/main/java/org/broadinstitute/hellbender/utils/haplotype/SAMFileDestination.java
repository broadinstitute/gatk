package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.SAMFileHeader;
import java.nio.file.Path;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;

/**
 * Class used to direct output from a HaplotypeBAMWriter to a bam/sam file.
 */
public final class SAMFileDestination extends HaplotypeBAMDestination {
    private final SAMFileGATKReadWriter samWriter;

    /**
     * Create a new file destination.
     *
     * @param outPath path where output is written (doesn't have to be local)
     * @param createBamOutIndex true to create an index file for the bamout
     * @param createBamOutMD5 true to create an md5 file for the bamout
     * @param sourceHeader SAMFileHeader used to seed the output SAMFileHeader for this destination, must not be null
     * @param haplotypeReadGroupID  read group ID used when writing haplotypes as reads
     */
    public SAMFileDestination(
            final Path outPath,
            final boolean createBamOutIndex,
            final boolean createBamOutMD5,
            final SAMFileHeader sourceHeader,
            final String haplotypeReadGroupID)
    {
        super(sourceHeader, haplotypeReadGroupID);
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
        samWriter.addRead(read);
    }
}