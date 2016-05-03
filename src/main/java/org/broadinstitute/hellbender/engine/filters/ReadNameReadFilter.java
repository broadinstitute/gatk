package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep only reads with this read name.
 * Matching is done by case-sensitive exact match.
 */
public final class ReadNameReadFilter implements ReadFilter, CommandLineFilter, Serializable {
    private static final long serialVersionUID = 1L;

    private final static String readNameArgName = "readName";
    @Argument(fullName = readNameArgName, shortName = "rn", doc="Keep only reads with this read name", optional=true)
    public String readName;

    @Override
    public boolean test( final GATKRead read ) {
        return read.getName() != null && read.getName().equals(readName);
    }

    @Override
    public String validate() {
        String message = null;
        if (readName == null || readName.length() < 1) {
            message = "requires a value for \"" + readNameArgName + "\"";
        }

        return message;
    }

}
