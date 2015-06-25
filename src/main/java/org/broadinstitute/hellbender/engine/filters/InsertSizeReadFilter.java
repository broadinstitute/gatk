package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;

/**
 * Keep reads that are within a given max insert size.
 */
public final class InsertSizeReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "maxInsertSize", shortName = "maxInsert", doc="Keep only read pairs with insert size at most equal to the given value", optional=true)
    public int maxInsertSize = 1000000;

    @Override
    public boolean test(final SAMRecord read) {
        if (!read.getReadPairedFlag()) {
            return true;
        }
        //Note insert size is negative if mate maps to lower position than read so we take absolute value.
        return Math.abs(read.getInferredInsertSize()) <= maxInsertSize;
    }
}
