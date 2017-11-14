package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SerializableIntFunction;
import org.broadinstitute.hellbender.utils.SerializableToIntFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.ToIntFunction;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * Created by valentin on 10/3/17.
 */
@DefaultSerializer(SVIntervalLocator.Serializer.class)
public final class SVIntervalLocator {

    private static final long serialVersionUID = 1L;

    private final List<String> contigByIndex;
    private final Map<String, Integer> indexByContig;



    private SVIntervalLocator(final List<String> contigByIndex) {
        this.contigByIndex = contigByIndex;
        final int numberOfContigs = contigByIndex.size();
        indexByContig = new HashMap<>(numberOfContigs);
        for (int i = 0; i < numberOfContigs; i++) {
            indexByContig.put(contigByIndex.get(i), i);
        }
    }

    /**
     * Creates a locator based on a {@link SAMSequenceDictionary}.
     * <p>
     *  This locator will resolve the contig index using {@link SAMSequenceDictionary#getSequenceIndex(String)}.
     * </p>
     */
    public static SVIntervalLocator of(final SAMSequenceDictionary dictionary) {
        return new SVIntervalLocator(Utils.nonNull(dictionary, "the input dictionary cannot be null")
                .getSequences().stream()
                .map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList()));
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
        final Integer index = indexByContig.get(loc.getContig());
        Utils.nonNull(index, "the input location has an unknown contig: "+ loc.getContig());
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
        final int index = svInterval.getContig();
        Utils.validateArg(index >= 0 && index < contigByIndex.size(), "the input sv-interval index points to unexistent index");
        return new SimpleInterval(contigByIndex.get(svInterval.getContig()), svInterval.getStart(), svInterval.getEnd() - 1);
    }

    public static class Serializer extends com.esotericsoftware.kryo.Serializer<SVIntervalLocator> {

        @Override
        public void write(final Kryo kryo, final Output output, final SVIntervalLocator object) {
            final int numberOfContigs = object.contigByIndex.size();
            output.writeInt(numberOfContigs);
            for (int i = 0; i < numberOfContigs; i++) {
                output.writeString(object.contigByIndex.get(i));
            }
        }

        @Override
        public SVIntervalLocator read(Kryo kryo, Input input, Class<SVIntervalLocator> type) {
            final int numberOfContigs = input.readInt();
            final List<String> contigsByIndex = new ArrayList<>(numberOfContigs);
            for (int i = 0; i < numberOfContigs; i++) {
                contigsByIndex.add(input.readString());
            }
            return new SVIntervalLocator(contigsByIndex);
        }
    }
}
