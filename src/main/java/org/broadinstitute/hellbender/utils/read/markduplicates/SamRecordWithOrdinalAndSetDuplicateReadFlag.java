package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SamRecordWithOrdinal;

/**
 * This class sets the duplicate read flag as the result state when examining sets of records.
 *
 * @author nhomer
 */
public final class SamRecordWithOrdinalAndSetDuplicateReadFlag extends SamRecordWithOrdinal {

    public SamRecordWithOrdinalAndSetDuplicateReadFlag() {
        super();
    }

    public SamRecordWithOrdinalAndSetDuplicateReadFlag(final SAMRecord record, final long recordIndex) {
        super(record, recordIndex);
    }

    @Override
    public void setResultState(final boolean resultState) {
        this.getRecord().setDuplicateReadFlag(resultState);
    }
}
