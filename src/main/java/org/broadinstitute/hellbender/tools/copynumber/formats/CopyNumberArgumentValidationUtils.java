package org.broadinstitute.hellbender.tools.copynumber.formats;

import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.OptionalInt;
import java.util.stream.IntStream;

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

    /**
     * Validate that a list of locatables is sorted according to a sequence dictionary and contains no duplicates or overlaps.
     */
    public static <T extends Locatable> void validateIntervals(final List<T> intervals,
                                                               final SAMSequenceDictionary sequenceDictionary) {
        Utils.nonNull(intervals);
        Utils.nonNull(sequenceDictionary);
        if (!Ordering.from(IntervalUtils.getDictionaryOrderComparator(sequenceDictionary)).isStrictlyOrdered(intervals)) {
            throw new IllegalArgumentException("Records were not strictly sorted in dictionary order.");
        }
        final OptionalInt failureIndex = IntStream.range(1, intervals.size())
                .filter(i -> IntervalUtils.overlaps(intervals.get(i - 1), intervals.get(i)))
                .findFirst();
        if (failureIndex.isPresent()) {
            final int index = failureIndex.getAsInt();
            throw new IllegalArgumentException(
                    String.format("Records contain at least two overlapping intervals: %s and %s",
                            intervals.get(index - 1), intervals.get(index)));
        }
    }
}
