package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;

/**
 * Keep only reads whose length is >= min value and <= max value.
 */
public final class ReadLengthReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "maxReadLength", shortName = "maxRead", doc="Keep only reads with length at most equal to the specified value", optional=false)
    public int maxReadLength;

    @Argument(fullName = "minReadLength", shortName = "minRead", doc="Keep only reads with length at least equal to the specified value", optional=false)
    public int minReadLength = 1;

    @Override
    public boolean test(final SAMRecord read) {
        return read.getReadLength() >= minReadLength && read.getReadLength() <= maxReadLength;
    }
}
