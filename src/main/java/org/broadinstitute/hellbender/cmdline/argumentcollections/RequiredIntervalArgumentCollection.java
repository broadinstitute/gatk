package org.broadinstitute.hellbender.cmdline.argumentcollections;


import org.broadinstitute.hellbender.cmdline.Argument;

import java.util.ArrayList;
import java.util.List;


/**
 * An ArgumentCollection that requires one or more intervals be specified with -L at the command line
 */
public final class RequiredIntervalArgumentCollection extends IntervalArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "intervals", shortName = "L", doc = "One or more genomic intervals over which to operate", optional = false)
    final protected List<String> intervalStrings = new ArrayList<>();

    @Override
    protected List<String> getIntervalStrings() {
        return intervalStrings;
    }

    @Override
    protected void addToIntervalStrings(String newInterval) {
        intervalStrings.add(newInterval);
    }
}
