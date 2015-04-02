package org.broadinstitute.hellbender.tools.readersplitters;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Splits a reader based on a value from a read group.
 * @param <T> The type of the value from the read group.
 */
public abstract class ReadGroupSplitter<T> extends ReaderSplitter<T> {
    @Override
    public T getSplitBy(final GATKRead record, final SAMFileHeader header) {
        return getSplitByFunction().apply(ReadUtils.getSAMReadGroupRecord(record, header));
    }

    @Override
    public List<T> getSplitsBy(final SAMFileHeader header) {
        return header.getReadGroups().stream().map(getSplitByFunction()).collect(Collectors.toList());
    }

    /**
     * @return A function that selects a value from a read group record.
     */
    protected abstract Function<SAMReadGroupRecord, T> getSplitByFunction();
}
