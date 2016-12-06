package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep reads that are within a given max fragment length.
 */
public final class FragmentLengthReadFilter extends ReadFilter implements Serializable  {
    private static final long serialVersionUID = 1l;

    @Argument(fullName = "maxFragmentLength",
            shortName = "maxFragmentLength",
            doc = "Keep only read pairs with fragment length at most equal to the given value",
            optional = true)
    public int maxFragmentLength = 1000000;

    @Override
    public boolean test( final GATKRead read ) {
        if ( ! read.isPaired() ) {
            return true;
        }
        //Note fragment length is negative if mate maps to lower position than read so we take absolute value.
        return Math.abs(read.getFragmentLength()) <= maxFragmentLength;
    }
}
