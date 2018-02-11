package org.broadinstitute.hellbender.tools.copynumber.utils;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.Set;

/**
 * A simple wrapper around {@link OverlapDetector} to provide naive caching and ensure that overlap sets
 * only contain a single interval.
 */
public final class CachedOverlapDetector<T extends Locatable> {
    private final OverlapDetector<T> overlapDetector;
    private T cachedResult;

    public CachedOverlapDetector(final List<T> intervals) {
        Utils.nonEmpty(intervals);
        this.overlapDetector = OverlapDetector.create(intervals);
        cachedResult = intervals.get(0);
    }

    /**
     * TODO enforce the locatable to be of size 1
     *
     * We check the previously cached result first.  Assuming that queries will be made in sorted order,
     * this may slightly save on lookup time.
     * @return {@code null} if no interval overlaps {@code locatable}
     */
    public T getOverlap(final Locatable locatable) {
        if (IntervalUtils.overlaps(cachedResult, locatable)) {
            return cachedResult;
        }
        final Set<T> overlaps = overlapDetector.getOverlaps(locatable);
        if (overlaps.size() > 1) {
            //should not reach here since intervals are checked for overlaps;
            //performing a redundant check to protect against future code changes
            throw new GATKException.ShouldNeverReachHereException("Intervals should be non-overlapping, " +
                    "so at most one interval should intersect with the center of a fragment.");
        }
        final T firstOverlap = overlaps.stream().findFirst().orElse(null);
        if (firstOverlap != null) {
            cachedResult = firstOverlap;
        }
        return firstOverlap;
    }

    public OverlapDetector<T> getOverlapDetector() {
        return overlapDetector;
    }
}
