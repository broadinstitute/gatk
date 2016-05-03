package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep reads that are within a given max fragment length.
 */
public final class FragmentLengthReadFilter implements ReadFilter, CommandLineFilter, Serializable {
    private static final long serialVersionUID = 1L;

    private final static String maxLengthArgName = "maxFragmentLength";

    @Argument(fullName = maxLengthArgName,
            shortName = "maxFragment",
            doc="Keep only read pairs with fragment length at most equal to the given value",
            optional=true)
    public int maxFragmentLength = 1000000;

    @Override
    public boolean test( final GATKRead read ) {
        if ( ! read.isPaired() ) {
            return true;
        }
        //Note fragment length is negative if mate maps to lower position than read so we take absolute value.
        return Math.abs(read.getFragmentLength()) <= maxFragmentLength;
    }

    @Override
    public String validate() {
        String message = null;
        if (maxFragmentLength < 1) {
            message = "requires a non-zero value for \"" + maxLengthArgName + "\"";
        }

        return message;
    }

}
