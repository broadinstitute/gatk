package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Predicate;
import java.util.stream.IntStream;

/**
 * Wrapper/adapter for {@link VariantFilter} that counts the number of variants filtered,
 * and provides a filter count summary.
 *
 * For filters with complex predicates that are composed of compound and/or operators, provides counts
 * by level for each component predicate. All counts reflect short-circuited evaluation of and/or operators
 * (ie., not all variant filters in a compound predicate will necessarily be evaluated every time). Also note
 * that the count for a compound predicate does not always equal the sum of the counts of it's component
 * predicates, i.e. an "or" filter will report a filter count of 1 in the case where both component predicates
 * are evaluated and both fail, but the the individual component filters will report a count of 1 at the next
 * level.
 * 
 */
public class CountingVariantFilter implements VariantFilter {

    private static final long serialVersionUID = 1L;

    // Underlying VariantFilter we delegate to if we're wrapping a simple VariantFilter.
    private final VariantFilter delegateFilter;

    // Number of variants filtered by this filter
    protected long filteredCount = 0;

    public CountingVariantFilter(final VariantFilter variantFilter) {
        Utils.nonNull(variantFilter);
        delegateFilter = variantFilter;
    }

    // Used only by the nested CountingBinopVariantFilter subclass and its derivatives, which must
    // override the test method with an implementation that does not depend on delegateFilter.
    private CountingVariantFilter() {
        delegateFilter = null;
    }

    // Return the number of variants filtered by this filter
    public long getFilteredCount() {
        return filteredCount;
    }

    public void resetFilteredCount() {
        filteredCount = 0;
    }

    public String getName() {return delegateFilter.getClass().getSimpleName();}

    // Returns a summary line with filter counts organized by level
    public String getSummaryLine() {return getSummaryLineForLevel(0);}

    protected String getSummaryLineForLevel(final int indentLevel) {
        if (0 == filteredCount) {
            return 0 == indentLevel ? "" : "No variants filtered by: " + getName();
        }
        else {
            return getIndentString(indentLevel) + Long.toString(filteredCount) + " variants(s) filtered by: " + getName() + " \n";
        }
    }

    protected String getIndentString(final int indentLevel) {
        final StringBuilder bldr = new StringBuilder();
        IntStream.range(0, indentLevel).forEach(i -> bldr.append("  "));
        return bldr.toString();
    }

    /**
     * Specialization of {@link #and(Predicate)} so that CountingVariantFilter and'ed with other CountingVariantFilter produce a CountingVariantFilter
     */
    //@Override
    public CountingVariantFilter and(final CountingVariantFilter other) {
        Utils.nonNull(other);
        return new CountingVariantFilter.CountingAndVariantFilter(this, other);
    }

    /**
     * Specialization of {@link #or(Predicate)} so that CountingVariantFilter ored with other CountingVariantFilter produce a CountingVariantFilter
     */
    //@Override
    public CountingVariantFilter or(final CountingVariantFilter other) {
        Utils.nonNull(other);
        return new CountingVariantFilter.CountingOrVariantFilter(this, other);
    }

    /**
     * Specialization of negate so that the resulting object is still a CountingVariantFilter
     */
    @Override
    public CountingVariantFilter negate() {
        return new CountingVariantFilter.CountingNegateVariantFilter(this);
    }

    @Override
    public boolean test(final VariantContext variant) {
        final boolean accept = delegateFilter.test(variant);
        if (!accept) {
            filteredCount++;
        }
        return accept;
    }

    private class CountingNegateVariantFilter extends CountingVariantFilter {
        private static final long serialVersionUID = 1L;

        CountingVariantFilter delegateCountingFilter;

        public CountingNegateVariantFilter(CountingVariantFilter delegate) {
            this.delegateCountingFilter = delegate;
        }

        @Override
        public boolean test(VariantContext variant) {
            final boolean accept = !delegateCountingFilter.test(variant);
            if (!accept) {
                filteredCount++;
            }
            return accept;
        }

        @Override
        public String getName() {
            return "Not " + delegateCountingFilter.getName();
        }
    }

    /**
     * Private class for Counting binary operator (and/or) filters; these keep track of how many variants are filtered at
     * each level of filter nesting.
     *
     * Subclasses must override the test method.
     */
    private abstract class CountingBinopVariantFilter extends CountingVariantFilter {

        private static final long serialVersionUID = 1L;

        final CountingVariantFilter lhs;
        final CountingVariantFilter rhs;

        public CountingBinopVariantFilter(final CountingVariantFilter lhs, final CountingVariantFilter rhs) {
            Utils.nonNull(lhs);
            Utils.nonNull(rhs);
            this.lhs = lhs;
            this.rhs = rhs;
        }

        @Override
        protected String getSummaryLineForLevel(final int indentLevel) {
            final String indent = getIndentString(indentLevel);
            if (0 == filteredCount) {
                return "No variants filtered by: " + getName();
            }
            else {
                return indent + Long.toString(filteredCount) + " variants(s) filtered by: " + getName() + "\n"
                        + (lhs.getFilteredCount() > 0 ? indent + lhs.getSummaryLineForLevel(indentLevel + 1) : "")
                        + (rhs.getFilteredCount() > 0 ? indent + rhs.getSummaryLineForLevel(indentLevel + 1) : "");
            }
        }

        @Override
        public abstract String getName();
    }

    /**
     * Private class for Counting AND filters
     */
    private final class CountingAndVariantFilter extends CountingVariantFilter.CountingBinopVariantFilter {

        private static final long serialVersionUID = 1L;

        private CountingAndVariantFilter(final CountingVariantFilter lhs, final CountingVariantFilter rhs) {
            super(lhs, rhs);
        }

        @Override
        public boolean test(final VariantContext variant) {
            final boolean accept = lhs.test(variant) && rhs.test(variant);
            if (!accept) {
                filteredCount++;
            }
            return accept;
        }

        @Override
        public String getName() {
            return "(" + lhs.getName() + " AND " + rhs.getName() + ")";
        }
    }

    /**
     * Private class for Counting OR filters
     */
    private final class CountingOrVariantFilter extends CountingVariantFilter.CountingBinopVariantFilter {

        private static final long serialVersionUID = 1L;

        private CountingOrVariantFilter(final CountingVariantFilter lhs, final CountingVariantFilter rhs) {
            super(lhs, rhs);
        }

        @Override
        public boolean test(final VariantContext variant) {
            final  boolean accept = lhs.test(variant) || rhs.test(variant);
            if (!accept) {
                filteredCount++;
            }
            return accept;
        }

        @Override
        public String getName() {
            return "(" + lhs.getName() + " OR " + rhs.getName() + ")";
        }
    }

}
