package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.SerializableFunction;

import java.util.Objects;
import java.util.function.UnaryOperator;

/**
 * Classes which perform transformations from SamRecord -> SAMRecord should implement this interface by overriding {@link SerializableFunction<SAMRecord,SAMRecord>#apply(SAMRecord)}
 */
@FunctionalInterface
public interface SamRecordTransformer extends UnaryOperator<SAMRecord>, SerializableFunction<SAMRecord, SAMRecord> {
    public static final long serialVersionUID = 1L;

    //HACK: These methods are a hack to get to get the type system to accept compositions of SAMRecords.
    @SuppressWarnings("overloads")
    default SamRecordTransformer andThen(final SamRecordTransformer after) {
        Objects.requireNonNull(after);
        return (SAMRecord r) -> after.apply(apply(r));
    }

    @SuppressWarnings("overloads")
    default SamRecordTransformer compose(final SamRecordTransformer before) {
        Objects.requireNonNull(before);
        return (SAMRecord r) -> apply(before.apply(r));
    }

    static SamRecordTransformer identity(){
        return read -> read;
    }
}
