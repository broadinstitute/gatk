package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;

/**
 * Keep only reads with this read name.
 * Matching is done by case-sensitive exact match.
 */
public final class ReadNameReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;

     @Argument(fullName = "readName", shortName = "rn", doc="Keep only reads with this read name", optional=false)
    public String readName;

    @Override
    public boolean test(final SAMRecord read) {
        return read.getReadName().equals(readName);
    }
}
