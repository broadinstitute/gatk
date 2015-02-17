package org.broadinstitute.hellbender.tools.walkers.bqsr;

import htsjdk.samtools.SAMRecord;

import java.util.function.Function;


public interface ReadTransformer extends Function<SAMRecord ,SAMRecord> {
    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadTransformers.
    default public ReadTransformer andThen(ReadTransformer after) {
        return Function.super.andThen(after)::apply;
    }

    default public  ReadTransformer compose(ReadTransformer before) {
        return Function.super.compose(before)::apply;
    }

    static public ReadTransformer identity(){
        return read -> read;
    }

}
