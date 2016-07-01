package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Only use reads from the specified read group.
 * Matching is done by case-sensitive exact match.
 */
public final class ReadGroupReadFilter implements ReadFilter, CommandLineFilter, Serializable {
    private static final long serialVersionUID = 1L;

    private final static String groupArgName = "readGroup";
    @Argument(fullName = groupArgName, shortName = "goodRG", doc="The name of the read group to keep", optional=true)
    public String readGroup = null;

    @Override
    public boolean test( final GATKRead read ) {
        final String rg = read.getReadGroup();
        return readGroup != null && rg.equals(this.readGroup);
    }

    @Override
    public String validate() {
        String message = null;
        if (readGroup == null || readGroup.length() < 1) {
            message = "requires a value for \"" + groupArgName + "\"";
        }

        return message;
    }

}