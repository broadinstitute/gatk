package org.broadinstitute.hellbender.tools.readersplitters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Splits a reader by some value.
 * @param <T> Type of the value that will be split by.
 */
public abstract class ReaderSplitter<T> {
    /**
     * Returns the values from the header that will be used to split the reader.
     * @param header The header of the reader.
     * @return The list of possibly values from the header for this splitter.
     */
    public abstract List<T> getSplitsBy(final SAMFileHeader header);

    /**
     * Returns the value from this record for this splitter.
     * @param record The record.
     * @param header Header for the record
     * @return The value from the record for this splitter.
     */
    public abstract T getSplitBy(final GATKRead record, final SAMFileHeader header);
}
