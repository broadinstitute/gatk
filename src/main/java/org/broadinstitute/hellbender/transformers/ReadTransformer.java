package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Objects;
import java.util.function.UnaryOperator;

/**
 * Classes which perform transformations from GATKRead -> GATKRead should implement this interface by overriding {@link SerializableFunction<GATKRead,GATKRead>#apply(GATKRead)}
 */
@FunctionalInterface
public interface ReadTransformer extends UnaryOperator<GATKRead>, SerializableFunction<GATKRead, GATKRead> {
    public static final long serialVersionUID = 1L;

    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadTransformers.
    @SuppressWarnings("overloads")
    default ReadTransformer andThen(ReadTransformer after) {
        Objects.requireNonNull(after);
        return (GATKRead r) -> after.apply(apply(r));
    }

    @SuppressWarnings("overloads")
    default ReadTransformer compose(ReadTransformer before) {
        Objects.requireNonNull(before);
        return (GATKRead r) -> apply(before.apply(r));
    }

    static ReadTransformer identity(){
        return read -> read;
    }
}
