package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep only reads whose length is >= min value and <= max value.
 */
public final class ReadLengthReadFilter implements ReadFilter, CommandLineFilter, Serializable {
    private static final long serialVersionUID = 1L;

    private static final String maxLengthArgName = "maxReadLength";

    @Argument(fullName = maxLengthArgName, shortName = "maxRead", doc="Keep only reads with length at most equal to the specified value", optional=false)
    public int maxReadLength = 0;

    private static final String minLengthArg = "minReadLength";
    @Argument(fullName = minLengthArg, shortName = "minRead", doc="Keep only reads with length at least equal to the specified value", optional=false)
    public int minReadLength = 0;

    @Override
    public boolean test( final GATKRead read ) {
        return read.getLength() >= minReadLength && read.getLength() <= maxReadLength;
    }

    @Override
    public String validate() {
        String message = null;
        if (maxReadLength <= 0 || minReadLength <= 0) {
            message = "requires non-zero argument values for \"" +
                    minLengthArg + "\" and \"" + maxLengthArgName + "\"";
        }

        return message;
    }

}
