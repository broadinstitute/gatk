package org.broadinstitute.hellbender.engine.filters;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import htsjdk.samtools.SAMRecord;

import java.io.Serializable;
import java.util.Objects;
import java.util.function.Function;
import java.util.function.Predicate;

@FunctionalInterface
public interface ReadFilter extends Predicate<SAMRecord>, SerializableFunction<SAMRecord, Boolean>{

    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadFilters.
    default ReadFilter and(ReadFilter other ) {
            Objects.requireNonNull(other);
            return (t) -> test(t) && other.test(t);
    }

    default ReadFilter or(ReadFilter other ) {
        Objects.requireNonNull(other);
        return (t) -> test(t) || other.test(t);
    }

    @Override
    default ReadFilter negate(){
        return (t) -> !test(t);
    }

    default Boolean apply(SAMRecord read){
        return test(read);
    }
}
