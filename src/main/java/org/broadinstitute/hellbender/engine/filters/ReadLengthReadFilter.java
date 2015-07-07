package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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
    public boolean test( final GATKRead read ) {
        return read.getLength() >= minReadLength && read.getLength() <= maxReadLength;
    }
}
