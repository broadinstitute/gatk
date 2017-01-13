package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

/**
 * Wrapper/adapter for {@link VariantFilter} that counts the number of variants filtered,
 * and provides a filter count summary.
 * <p>
 * For filters with complex predicates that are composed of compound and/or operators, provides counts
 * by level for each component predicate. All counts reflect short-circuited evaluation of and/or operators
 * (ie., not all variant filters in a compound predicate will necessarily be evaluated every time). Also note
 * that the count for a compound predicate does not always equal the sum of the counts of it's component
 * predicates, i.e. an "or" filter will report a filter count of 1 in the case where both component predicates
 * are evaluated and both fail, but the the individual component filters will report a count of 1 at the next
 * level.
 */
public class CountingVariantFilter implements VariantFilter {
    private static final long serialVersionUID = 1L;

    private static final String objectName = "variant";

    private CountingFilter<VariantContext> countingFilter;

    public CountingVariantFilter(final VariantFilter variantFilter) {
        this(new CountingFilter<>(variantFilter, objectName));
    }

    private CountingVariantFilter(final CountingFilter<VariantContext> countingFilter) {
        this.countingFilter = countingFilter;
    }

    /**
     * Return a composite (and) {@code CountingVariantFilter} constructed from a list of
     * {@link org.broadinstitute.hellbender.engine.filters.VariantFilter}.
     * The resulting filter honors the order of the input list
     * and tests the filter conditions in the same order as the iteration order of the input list.
     * @param variantFilters If null or empty, the ALLOW_ALL_VARIANTS variant filter will be returned
     * @return Composite CountingVariantFilter
     */
    public static CountingVariantFilter fromList(final List<VariantFilter> variantFilters) {
        if (variantFilters == null || variantFilters.isEmpty()) {
            return new CountingVariantFilter(VariantFilterLibrary.ALLOW_ALL_VARIANTS);
        }
        CountingVariantFilter compositeFilter = new CountingVariantFilter(variantFilters.get(0));
        for (int i = 1; i < variantFilters.size(); i++) {
            compositeFilter = compositeFilter.and(new CountingVariantFilter(variantFilters.get(i)));
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

    @Override
    public boolean test(final VariantContext variantContext) {
        return countingFilter.test(variantContext);
    }

    @Override
    public VariantFilter and(VariantFilter filter) {
        return new CountingVariantFilter(this.countingFilter.and(new CountingFilter<>(filter, objectName)));
    }

    @Override
    public VariantFilter or(VariantFilter filter) {
        return new CountingVariantFilter(this.countingFilter.or(new CountingFilter<>(filter, objectName)));
    }

    @Override
    public CountingVariantFilter negate() {
        return new CountingVariantFilter(this.countingFilter.negate());
    }

    public CountingVariantFilter and(CountingVariantFilter filter) {
        return new CountingVariantFilter(this.countingFilter.and(filter.countingFilter));
    }

    public CountingVariantFilter or(CountingVariantFilter filter) {
        return new CountingVariantFilter(this.countingFilter.or(filter.countingFilter));
    }
}
