package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
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

    /**
     * Return a composite (and) {@code ReadTransformer} constructed from a list of {@code ReadTransformer}. Each
     * transformer in the list is first initialized with the {@code SAMFileHeader} param. The resulting transformer
     * honors the order of the input list and tests the transformer conditions in the same order as the iteration
     * order of the input list.
     * @param readTransformers If null or empty, the ALLOW_ALL_READS read transformer will be returned
     * @param samHeader {@code SAMFileHeader} used to initialize each transformer. May not be null
     * @return Composite ReadTransformer
     */
    static ReadTransformer fromList(final List<ReadTransformer> readTransformers, final SAMFileHeader samHeader) {
        Utils.nonNull(samHeader, "SAMFileHeader must not be null");
        if (readTransformers == null || readTransformers.isEmpty()) {
            return ReadTransformer.identity();
        }
        ReadTransformer compositeTransformer = readTransformers.get(0);
        for (int i = 1; i < readTransformers.size(); i++) {
            compositeTransformer = compositeTransformer.andThen(readTransformers.get(i));
        }
        return compositeTransformer;
    }



}
