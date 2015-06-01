package org.broadinstitute.hellbender.engine.filters;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import htsjdk.samtools.SAMRecord;

import java.io.Serializable;
import java.util.Objects;
import java.util.function.Function;
import java.util.function.Predicate;

@FunctionalInterface
/**
 * Filters which operate on {@link SAMRecord} should implement this interface by overriding {@link #test(SAMRecord)}
 *
 * ReadFilter extends Predicate and SerializableFunction.  It provides a default implementation of apply based on the
 * implmenting class's implementation of test().
 */
public interface ReadFilter extends Predicate<SAMRecord>, SerializableFunction<SAMRecord, Boolean>{

    // It turns out, this is necessary. Please don't remove it.
    // Without this line, we see the following error:
    // java.io.InvalidClassException: org.broadinstitute.hellbender.engine.filters.ReadFilter; local class incompatible:
    // stream classdesc serialVersionUID = -5040289903122017748, local class serialVersionUID = 6814309376393671214
    static final long serialVersionUID = 1L;

    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadFilters.
    /**
     * Specialization of {@link #and(Predicate)} so that ReadFilters anded with other ReadFilters produce a ReadFilter
     */
    default ReadFilter and(ReadFilter other ) {
            Objects.requireNonNull(other);
            return (t) -> test(t) && other.test(t);
    }

    /**
     * Specialization of {@link #or(Predicate)} so that ReadFilters ored with other ReadFilters produce a ReadFilter
     */
    default ReadFilter or(ReadFilter other ) {
        Objects.requireNonNull(other);
        return (t) -> test(t) || other.test(t);
    }

    /**
     * Specialization of negate so that the resulting object is still a ReadFilter
     */
    @Override
    default ReadFilter negate(){
        return (t) -> !test(t);
    }

    @Override
    default Boolean apply(SAMRecord read){
        return test(read);
    }

    @Override
    boolean test(SAMRecord read);
}
