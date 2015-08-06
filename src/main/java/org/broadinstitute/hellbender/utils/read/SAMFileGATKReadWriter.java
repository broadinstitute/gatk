package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileWriter;

/**
 * A GATKRead writer that writes to a SAM/BAM file.
 *
 * Converts each read to SAMRecord in the process, which may be a lossy operation if the
 * read is not already in SAM format.
 */
public final class SAMFileGATKReadWriter implements GATKReadWriter {

    private final SAMFileWriter samWriter;

    public SAMFileGATKReadWriter( final SAMFileWriter samWriter ) {
        this.samWriter = samWriter;
    }

    @Override
    public void addRead( GATKRead read ) {
        samWriter.addAlignment(read.convertToSAMRecord(samWriter.getFileHeader()));
    }

    @Override
    public void close() {
        samWriter.close();
    }
}
