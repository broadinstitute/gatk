package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Only use reads from the specified read group.
 * Matching is done by case-sensitive exact match.
 */
public final class ReadGroupReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "read_group_to_keep", shortName = "goodRG", doc="The name of the read group to keep", optional=false)
    public String readGroup = null;

    @Override
    public boolean test( final GATKRead read ) {
        final String rg = read.getReadGroup();
        return readGroup != null && rg.equals(this.readGroup);
    }
}