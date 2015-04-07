package org.broadinstitute.hellbender.engine.filters;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import htsjdk.samtools.SAMRecord;

import java.io.Serializable;
import java.util.function.Function;
import java.util.function.Predicate;

@FunctionalInterface
public interface ReadFilter extends Predicate<SAMRecord>, SerializableFunction<SAMRecord, Boolean>, Serializable{

    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadFilters.
    default ReadFilter and(ReadFilter filter ) {
        return Predicate.super.and(filter)::test;
    }

    default ReadFilter or(ReadFilter filter ) {
        return Predicate.super.or(filter)::test;
    }

    default ReadFilter negate(){
        return Predicate.super.negate()::test;
    }

    default Boolean apply(SAMRecord read){
        return test(read);
    }
}
