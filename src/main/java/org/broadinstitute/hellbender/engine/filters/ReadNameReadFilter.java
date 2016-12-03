package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep only reads with this read name.
 * Matching is done by case-sensitive exact match.
 */
public final class ReadNameReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "readName", shortName = "readName", doc="Keep only reads with this read name", optional=false)
    public String readName = null;

    @Override
    public boolean test( final GATKRead read ) {
        return read.getName() != null && read.getName().equals(readName);
    }

}
