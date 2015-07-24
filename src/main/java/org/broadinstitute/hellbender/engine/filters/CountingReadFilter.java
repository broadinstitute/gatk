package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Predicate;
import java.util.stream.IntStream;

/**
 * Wrapper/adapter for ReadFilter that counts the number of reads filtered, and provides a filter count summary.
 *
 * For filters with complex predicates that are composed of compound and/or operators, provides counts
 * by level for each component predicate. All counts reflect short-circuited evaluation of and/or operators
 * (ie., not all read filters in a compound predicate will necessarily be evaluated every time). Also note
 * that the count for a compound predicate does not always equal the sum of the counts of it's component
 * predicates, i.e. an "or" filter will report a filter count of 1 in the case where both component predicates
 * are evaluated and both fail, but the the individual component filters will report a count of 1 at the next
 * level.
 *
 */
public class CountingReadFilter implements ReadFilter {

    private static final long serialVersionUID = 1L;

    // Underlying ReadFilter we delegate to if we're wrapping a simple ReadFilter.
    private final ReadFilter delegateFilter;

    // Name for this filter, used to report filter count summaries
    protected final String name;

    // Number of reads filtered by this filter
    protected long filteredCount = 0;

    public CountingReadFilter(final String filterName, final ReadFilter readFilter) {
        Utils.nonNull(readFilter);
        name = filterName;
        delegateFilter = readFilter;
    }

    // Used only by the nested CountingBinopReadFilter subclass and its derivatives, which must
    // override the test method with an implementation that does not depend on delegateFilter.
    private CountingReadFilter(final String filterName) {
        name = filterName;
        delegateFilter = null;
    }

    // Return the number of reads filtered by this filter
    public long getFilteredCount() {
        return filteredCount;
    }

    public void resetFilteredCount() {
        filteredCount = 0;
    }

    public String getName() {return name;}

    // Returns a summary line with filter counts organized by level
    public String getSummaryLine() {return getSummaryLineForLevel(0);}

    protected String getSummaryLineForLevel(final int indentLevel) {
        if (0 == filteredCount) {
            return 0 == indentLevel ? "" : "No reads filtered by: " + getName();
        }
        else {
            return getIndentString(indentLevel) + Long.toString(filteredCount) + " read(s) filtered by: " + getName() + " \n";
        }
    }

    protected String getIndentString(final int indentLevel) {
        final StringBuilder bldr = new StringBuilder();
        IntStream.range(0, indentLevel).forEach(i -> bldr.append("  "));
        return bldr.toString();
    }

    /**
     * Specialization of {@link #and(Predicate)} so that CountingReadFilter and'ed with other CountingReadFilter produce a CountingReadFilter
     */
    //@Override
    public CountingReadFilter and(final CountingReadFilter other) {
        Utils.nonNull(other);
        return new CountingAndReadFilter(this, other);
    }

    /**
     * Specialization of {@link #or(Predicate)} so that CountingReadFilter ored with other CountingReadFilter produce a CountingReadFilter
     */
    //@Override
    public CountingReadFilter or(final CountingReadFilter other) {
        Utils.nonNull(other);
        return new CountingOrReadFilter(this, other);
    }

    /**
     * Specialization of negate so that the resulting object is still a CountingReadFilter
     */
    @Override
    public CountingReadFilter negate() {
        return new CountingReadFilter("(Not '" + getName() + "')",(t) -> !test(t));
    }

    @Override
    public boolean test(final GATKRead read) {
        final boolean accept = delegateFilter.test(read);
        if (!accept) {
            filteredCount++;
        }
        return accept;
    }

    /**
     * Private class for Counting binary operator (and/or) filters; these keep track of how many reads are filtered at
     * each level of filter nesting.
     *
     * Subclasses must override the test method.
     */
    private abstract class CountingBinopReadFilter extends CountingReadFilter {

        private static final long serialVersionUID = 1L;

        final CountingReadFilter lhs;
        final CountingReadFilter rhs;

        public CountingBinopReadFilter(final String name, final CountingReadFilter lhs, final CountingReadFilter rhs) {
            super(name);
            Utils.nonNull(lhs);
            Utils.nonNull(rhs);
            this.lhs = lhs;
            this.rhs = rhs;
        }

        @Override
        protected String getSummaryLineForLevel(final int indentLevel) {
            final String indent = getIndentString(indentLevel);
            if (0 == filteredCount) {
                return "No reads filtered by: " + getName();
            }
            else {
                return indent + Long.toString(filteredCount) + " read(s) filtered by: " + getName() + "\n"
                        + (lhs.getFilteredCount() > 0 ? indent + lhs.getSummaryLineForLevel(indentLevel + 1) : "")
                        + (rhs.getFilteredCount() > 0 ? indent + rhs.getSummaryLineForLevel(indentLevel + 1) : "");
            }
        }
    }

    /**
     * Private class for Counting AND filters
     */
    private class CountingAndReadFilter extends CountingBinopReadFilter {

        private static final long serialVersionUID = 1L;

        private CountingAndReadFilter(final CountingReadFilter lhs, final CountingReadFilter rhs) {
            super("(" + lhs.getName() + " AND " + rhs.getName() + ")", lhs, rhs);
        }

        @Override
        public boolean test(final GATKRead read) {
            final boolean accept = lhs.test(read) && rhs.test(read);
            if (!accept) {
                filteredCount++;
            }
            return accept;
        }
    }

    /**
     * Private class for Counting OR filters
     */
    private class CountingOrReadFilter extends CountingBinopReadFilter {

        private static final long serialVersionUID = 1L;

        private CountingOrReadFilter(final CountingReadFilter lhs, final CountingReadFilter rhs) {
            super("(" + lhs.getName() + " OR " + rhs.getName() + ")", lhs, rhs);
        }

        @Override
        public boolean test(final GATKRead read) {
            final  boolean accept = lhs.test(read) || rhs.test(read);
            if (!accept) {
                filteredCount++;
            }
            return accept;
        }
    }
}
