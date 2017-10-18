package org.broadinstitute.hellbender.tools.copynumber.formats;

import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyNumberArgumentValidationUtils {
    private CopyNumberArgumentValidationUtils() {}

    /**
     * Validate that the interval-argument collection parameters minimally modify the input intervals.
     */
    public static void validateIntervalArgumentCollection(final IntervalArgumentCollection intervalArgumentCollection) {
        Utils.validateArg(intervalArgumentCollection.getIntervalSetRule() == IntervalSetRule.UNION,
                "Interval set rule must be set to UNION.");
        Utils.validateArg(intervalArgumentCollection.getIntervalExclusionPadding() == 0,
                "Interval exclusion padding must be set to 0.");
        Utils.validateArg(intervalArgumentCollection.getIntervalPadding() == 0,
                "Interval padding must be set to 0.");
        Utils.validateArg(intervalArgumentCollection.getIntervalMergingRule() == IntervalMergingRule.OVERLAPPING_ONLY,
                "Interval merging rule must be set to OVERLAPPING_ONLY.");
    }
}
