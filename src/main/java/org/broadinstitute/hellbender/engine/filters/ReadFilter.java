package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.function.Predicate;

/**
 * Filters which operate on {@link GATKRead} should subclass this by overriding {@link #test(GATKRead)}
 *
 * ReadFilter implements Predicate and Serializable.  It provides a default implementation of apply
 * based on the implementing class's implementation of test().
 *
 * To be accessible from the command line, subclasses must have a zero-arg constructor and contain
 * (optional) ArgumentCollection or Argument annotated fields for filter-specific arguments. Arguments
 * that are optional=true should have an initial value since the command line parser will not require
 * the user to provide one; arguments that are optional=false should not have an initial value (and
 * used boxed types with null default/initial values) since otherwise they will be seen by the command
 * line parser as having been provided, and the user will not be notified that the value must be provided
 * on the command line).
 */
public abstract class ReadFilter implements Predicate<GATKRead>, Serializable {

    protected SAMFileHeader samHeader = null;

    protected final ReadFilter delegate;

    protected ReadFilter() { delegate = null; };
    protected ReadFilter(ReadFilter delegate) { this.delegate = delegate; }

    public void setHeader(SAMFileHeader samHeader) { this.samHeader = samHeader; }

    private class ReadFilterNegate extends ReadFilter {
        private static final long serialVersionUID = 1L;

        protected ReadFilterNegate(ReadFilter delegate) {super(delegate);}

        @Override
        public boolean test( GATKRead read ) {
            return !delegate.test(read);
        }
    }

    protected static class ReadFilterAnd extends ReadFilter {
        private static final long serialVersionUID = 1L;

        //package accesible for testing
        final protected ReadFilter other;

        public ReadFilterAnd(ReadFilter lhs, ReadFilter rhs) { super(lhs); this.other = rhs;}

        @Override
        public boolean test( GATKRead read ) { return delegate.test(read) && other.test(read);}
    }

    protected static class ReadFilterOr extends ReadFilter {
        private static final long serialVersionUID = 1L;

        final private ReadFilter other;

        public ReadFilterOr(ReadFilter lhs, ReadFilter rhs) { super(lhs); this.other = rhs;}

        @Override
        public boolean test( GATKRead read ) { return delegate.test(read) || other.test(read);}
    }

    // It turns out, this is necessary. Please don't remove it.
    // Without this line, we see the following error:
    // java.io.InvalidClassException: org.broadinstitute.hellbender.engine.filters.ReadFilter; local class incompatible:
    // stream classdesc serialVersionUID = -5040289903122017748, local class serialVersionUID = 6814309376393671214
    private static final long serialVersionUID = 1L;

    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadFilters.
    /**
     * Specialization of {@link #and(Predicate)} so that ReadFilters anded with other ReadFilters produce a ReadFilter
     */
    public ReadFilter and( ReadFilter other ) {
        Utils.nonNull(other);
        return new ReadFilterAnd(this, other);
    }

    /**
     * Specialization of {@link #or(Predicate)} so that ReadFilters ored with other ReadFilters produce a ReadFilter
     */
    public ReadFilter or( ReadFilter other ) {
        Utils.nonNull(other);
        return new ReadFilterOr(this, other);
    }

    /**
     * Specialization of negate so that the resulting object is still a ReadFilter
     */
    @Override
    public ReadFilter negate(){ return new ReadFilter.ReadFilterNegate(this); }

    @Override
    public abstract boolean test( GATKRead read );
}
