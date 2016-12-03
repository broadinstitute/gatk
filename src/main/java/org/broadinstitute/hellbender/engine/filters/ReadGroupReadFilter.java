package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Only use reads from the specified read group.
 * Matching is done by case-sensitive exact match.
 */
public final class ReadGroupReadFilter extends ReadFilter implements Serializable{
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "keepReadGroup", shortName = "keepReadGroup", doc="The name of the read group to keep", optional=false)
    public String readGroup = null;

    @Override
    public boolean test( final GATKRead read ) {
        final String rg = read.getReadGroup();
        return readGroup != null && rg.equals(this.readGroup);
    }
}