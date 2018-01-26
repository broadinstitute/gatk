
package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;

/**
 * An interval argument class that allows -L to be specified but does not require it.
 */
public final class OptionalIntervalArgumentCollection extends IntervalArgumentCollection {
    private static final long serialVersionUID = 1L;

    // Interval list files such as Picard interval lists are structured and require specialized parsing that
    // is handled by IntervalUtils, so use suppressFileExpansion to bypass command line parser auto-expansion.
    @Argument(fullName = StandardArgumentDefinitions.INTERVALS_LONG_NAME, shortName = StandardArgumentDefinitions.INTERVALS_SHORT_NAME,
            suppressFileExpansion = true, doc = "One or more genomic intervals over which to operate", optional = true)
    protected final List<String> intervalStrings = new ArrayList<>();

    @Override
    protected List<String> getIntervalStrings() {
        return intervalStrings;
    }

    @Override
    protected void addToIntervalStrings(String newInterval) {
        Utils.validate(traversalParameters == null, "addToIntervalStrings() cannot be called after interval parsing is complete");
        intervalStrings.add(newInterval);
    }
}

