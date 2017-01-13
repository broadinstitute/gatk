package org.broadinstitute.hellbender.engine.filters;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.function.Predicate;

/**
 * Wrapper/adapter for {@link org.broadinstitute.hellbender.engine.filters.ReadFilter} that counts the number of
 * reads filtered, and provides a filter count summary.
 *
 * For filters with complex predicates that are composed of compound and/or operators, provides counts
 * by level for each component predicate. All counts reflect short-circuited evaluation of and/or operators
 * (ie., not all read filters in a compound predicate will necessarily be evaluated every time). Also note
 * that the count for a compound predicate does not always equal the sum of the counts of it's component
 * predicates, i.e. an "or" filter will report a filter count of 1 in the case where both component predicates
 * are evaluated and both fail, but the the individual component filters will report a count of 1 at the next
 * level.
 */
public class CountingReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;

    private static final String objectName = "read";

    // Underlying ReadFilter we delegate to if we're wrapping a simple ReadFilter.
    @VisibleForTesting
    protected final CountingFilter<GATKRead> countingFilter;

    public CountingReadFilter(final ReadFilter readFilter) {
        Utils.nonNull(readFilter);
        countingFilter = new CountingFilter<>(readFilter, objectName);
    }

    private CountingReadFilter(final CountingFilter<GATKRead> countingFilter) {
        this.countingFilter = countingFilter;
    }

    /**
     * Return a composite (and) {@code CountingReadFilter} constructed from a list of
     * {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}. Each filter in the list is first
     * initialized with the {@code SAMFileHeader} param. The resulting filter honors the order of the input list
     * and tests the filter conditions in the same order as the iteration order of the input list.
     * @param readFilters If null or empty, the ALLOW_ALL_READS read filter will be returned
     * @param samHeader {@code SAMFileHeader} used to initialize each filter. May not be null
     * @return Composite CountingReadFilter
     */
    public static CountingReadFilter fromList(final List<ReadFilter> readFilters, final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader, "SAMFileHeader must not be null");
        if (readFilters == null || readFilters.isEmpty()) {
            return new CountingReadFilter(ReadFilterLibrary.ALLOW_ALL_READS);
        }
        readFilters.forEach(f -> f.setHeader(samHeader));
        CountingReadFilter compositeFilter = new CountingReadFilter(readFilters.get(0));
        for (int i = 1; i < readFilters.size(); i++) {
            compositeFilter = compositeFilter.and(new CountingReadFilter(readFilters.get(i)));
        }
        return compositeFilter;
    }

    // Return the number of reads filtered by this filter
    public long getFilteredCount() {
        return countingFilter.getFilteredCount();
    }

    public void resetFilteredCount() {
        countingFilter.resetFilteredCount();
    }

    public String getName() {return countingFilter.getClass().getSimpleName();}

    // Returns a summary line with filter counts organized by level
    public String getSummaryLine() {return countingFilter.getSummaryLine();}

    /**
     * Specialization of {@link #and(Predicate)} so that CountingReadFilter and'ed with other CountingReadFilter produce a CountingReadFilter
     */
    //@Override
    public CountingReadFilter and(final CountingReadFilter other) {
        Utils.nonNull(other);
        return new CountingReadFilter(this.countingFilter.and(new CountingFilter<>(other.countingFilter, objectName)));
    }

    /**
     * Specialization of {@link #or(Predicate)} so that CountingReadFilter ored with other CountingReadFilter produce a CountingReadFilter
     */
    //@Override
    public CountingReadFilter or(final CountingReadFilter other) {
        Utils.nonNull(other);
        return new CountingReadFilter(this.countingFilter.or(new CountingFilter<>(other.countingFilter, objectName)));
    }

    /**
     * Specialization of negate so that the resulting object is still a CountingReadFilter
     */
    @Override
    public CountingReadFilter negate() {
        return new CountingReadFilter(this.countingFilter.negate());
    }

    @Override
    public boolean test(final GATKRead read) {
        return countingFilter.test(read);
    }
}
