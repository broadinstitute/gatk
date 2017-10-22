package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SerializableIntFunction;
import org.broadinstitute.hellbender.utils.SerializableToIntFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.function.Function;
import java.util.function.ToIntFunction;
import java.util.stream.Collector;

/**
 * Created by valentin on 10/3/17.
 */
public final class SVIntervalLocator implements Serializable {

    private static final long serialVersionUID = 1L;

    private final SerializableToIntFunction<String> contigToIndex;
    private final SerializableIntFunction<String> indexToContig;

    public SVIntervalLocator(final SerializableToIntFunction<String> contigToIndex, final SerializableIntFunction<String> indexToContig) {
        this.contigToIndex = Utils.nonNull(contigToIndex);
        this.indexToContig = Utils.nonNull(indexToContig);
    }

    /**
     * Creates a locator based on a {@link SAMSequenceDictionary}.
     * <p>
     *  This locator will resolve the contig index using {@link SAMSequenceDictionary#getSequenceIndex(String)}.
     * </p>
     */
    public SVIntervalLocator(final SAMSequenceDictionary dictionary) {
        this(Utils.nonNull(dictionary, "the input dictionary cannot be null")::getSequenceIndex,
                idx -> dictionary.getSequence(idx).getSequenceName());
    }

    /**
     * Returns the {@link SVInterval} that represents the coordinates of a {@link Locatable}.
     * @param loc
     * @return never {@code null}.
     */
    public SVInterval toSVInterval(final Locatable loc) {
        return toSVInterval(loc, 0);
    }

    public SVInterval toSVInterval(final Locatable loc, final int padding) {
        Utils.nonNull(loc);
        ParamUtils.isPositiveOrZero(padding, "padding must 0 or positive");
        final int index = contigToIndex.applyAsInt(loc.getContig());
        final int start = Math.max(1, loc.getStart() - padding);
        final int end = loc.getEnd() + 1 + padding;
        return new SVInterval(index, start, end);
    }

    /**
     * Returns a function that returns the {@link SVInterval} that represent the coordinates of a {@link Locatable}.
     * @param <L>
     * @return
     */
    public <L extends Locatable> SerializableFunction<L, SVInterval> toSVIntervalLambda() {
        return this::toSVInterval;
    }

    public <L extends Locatable, V> Collector<L, ?, SVIntervalTree<V> > toTreeCollector(final Function<L, V> toValue) {
        return Collector.of(SVIntervalTree<V>::new,
                            (tree, loc) -> tree.put(toSVInterval(loc), toValue.apply(loc)),
                            (tree1, tree2) -> { tree1.putAll(tree2); return tree1; });
    }

    public SimpleInterval toSimpleInterval(final SVInterval svInterval) {
        Utils.nonNull(svInterval);
        return new SimpleInterval(indexToContig.apply(svInterval.getContig()), svInterval.getStart(), svInterval.getEnd() - 1);
    }
}
