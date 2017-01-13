package org.broadinstitute.hellbender.engine.filters;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Predicate;
import java.util.stream.IntStream;

/**
 * Wrapper/adapter for {@link Predicate<T>} that counts the number of objets filtered, and provides a filter count summary.
 *
 * For filters with complex predicates that are composed of compound and/or operators, provides counts
 * by level for each component predicate. All counts reflect short-circuited evaluation of and/or operators
 * (ie., not all filters in a compound predicate will necessarily be evaluated every time). Also note
 * that the count for a compound predicate does not always equal the sum of the counts of it's component
 * predicates, i.e. an "or" filter will report a filter count of 1 in the case where both component predicates
 * are evaluated and both fail, but the the individual component filters will report a count of 1 at the next
 * level.
 */
public class CountingFilter<T> implements Predicate<T> {
    private static final long serialVersionUID = 1L;

    // Underlying ReadFilter we delegate to if we're wrapping a simple ReadFilter.
    @VisibleForTesting
    protected final Predicate<T> delegateFilter;

    private final String objectName;

    // Number of objects filtered by this filter
    protected long filteredCount = 0;

    public CountingFilter(final Predicate<T> filter, final String objectName) {
        Utils.nonNull(filter);
        delegateFilter = filter;
        this.objectName = objectName;
    }

    // Used only by the nested CountingBinopReadFilter subclass and its derivatives, which must
    // override the test method with an implementation that does not depend on countingFilter.
    private CountingFilter() {
        delegateFilter = null;
        objectName = null;
    }

    // Return the number of objects filtered by this filter
    public long getFilteredCount() {
        return filteredCount;
    }

    public void resetFilteredCount() {
        filteredCount = 0;
    }

    public String getName() {return delegateFilter.getClass().getSimpleName();}

    protected String getObjectName() {return objectName;}

    // Returns a summary line with filter counts organized by level
    public String getSummaryLine() {return getSummaryLineForLevel(0);}

    protected String getSummaryLineForLevel(final int indentLevel) {
        if (0 == filteredCount) {
            return 0 == indentLevel ? "" : "No " + getObjectName() + " filtered by: " + getName();
        }
        else {
            return getIndentString(indentLevel) + Long.toString(filteredCount) + " " + getObjectName() + "(s) filtered by: " + getName() + " \n";
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
    // @Override
    public CountingFilter<T> and(final CountingFilter<T> other) {
        Utils.nonNull(other);
        return new CountingAndReadFilter<>(this, other);
    }

    /**
     * Specialization of {@link #or(Predicate)} so that CountingReadFilter ored with other CountingReadFilter produce a CountingReadFilter
     */
    //@Override
    public CountingFilter<T> or(final CountingFilter<T> other) {
        Utils.nonNull(other);
        return new CountingOrReadFilter<>(this, other);
    }

    /**
     * Specialization of negate so that the resulting object is still a CountingReadFilter
     */
    @Override
    public CountingFilter<T> negate() {
        return new CountingNegateReadFilter<>(this);
    }

    @Override
    public boolean test(final T t) {
        final boolean accept = delegateFilter.test(t);
        if (!accept) {
            filteredCount++;
        }
        return accept;
    }

    private static class CountingNegateReadFilter<T> extends CountingFilter<T> {
        private static final long serialVersionUID = 1L;

        CountingFilter<T> delegateCountingFilter;

        public CountingNegateReadFilter(CountingFilter<T> delegate) {
            this.delegateCountingFilter = delegate;
        }

        @Override
        protected final String getObjectName() {
            return delegateCountingFilter.getObjectName();
        }

        @Override
        public boolean test(T t) {
            final boolean accept = !delegateCountingFilter.test(t);
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
     * Private class for Counting binary operator (and/or) filters; these keep track of how many objects are filtered at
     * each level of filter nesting.
     *
     * Subclasses must override the test method.
     */
    private static abstract class CountingBinopReadFilter<T> extends CountingFilter<T> {

        private static final long serialVersionUID = 1L;

        protected final CountingFilter<T> lhs;
        protected final CountingFilter<T> rhs;

        public CountingBinopReadFilter(final CountingFilter<T> lhs, final CountingFilter<T> rhs) {
            Utils.nonNull(lhs);
            Utils.nonNull(rhs);
            Utils.validateArg((lhs.getObjectName().equals(rhs.getObjectName())), "different object names");
            this.lhs = lhs;
            this.rhs = rhs;
        }

        @Override
        protected final String getObjectName() {
            return lhs.getObjectName();
        }

        @Override
        protected String getSummaryLineForLevel(final int indentLevel) {
            final String indent = getIndentString(indentLevel);
            if (0 == filteredCount) {
                return "No " + lhs.getObjectName() + " filtered by: " + getName();
            }
            else {
                return indent + Long.toString(filteredCount) + " " + getObjectName() + "(s) filtered by: " + getName() + "\n"
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
    @VisibleForTesting
    protected static final class CountingAndReadFilter<T> extends CountingBinopReadFilter<T> {

        private static final long serialVersionUID = 1L;

        private CountingAndReadFilter(final CountingFilter<T> lhs, final CountingFilter<T> rhs) {
            super(lhs, rhs);
        }

        @Override
        public boolean test(final T t) {
            final boolean accept = lhs.test(t) && rhs.test(t);
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
    private static final class CountingOrReadFilter<T> extends CountingBinopReadFilter<T> {

        private static final long serialVersionUID = 1L;

        private CountingOrReadFilter(final CountingFilter<T> lhs, final CountingFilter<T> rhs) {
            super(lhs, rhs);
        }

        @Override
        public boolean test(final T t) {
            final  boolean accept = lhs.test(t) || rhs.test(t);
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
