package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.function.Function;
import java.util.stream.Collector;

/**
 * Class to translate {@link SVInterval} into {@link Locatable Locatables} and vice-versa
 */
@DefaultSerializer(SVIntervalLocator.Serializer.class)
public final class SVIntervalLocator {

    private final String[] contigByIndex;
    private final Object2IntMap<String> indexByContig;
    private final int[] contigLengths;

    @FunctionalInterface
    public interface LocatableFactory<L extends Locatable> {
        L create(final String contig, final int start, final int end);
    }

    private SVIntervalLocator(final String[] contigByIndex, final int[] contigLengths) {
        this.contigByIndex = contigByIndex;
        this.contigLengths = contigLengths;
        indexByContig = new Object2IntOpenHashMap<>(contigByIndex.length);
        for (int i = 0; i < contigByIndex.length; i++) {
            indexByContig.put(contigByIndex[i], i);
        }
    }

    /**
     * Creates a locator based on a {@link SAMSequenceDictionary}.
     * <p>
     *  This locator will resolve the contig index using
     *  {@link SAMSequenceDictionary#getSequenceIndex(String)}.
     * </p>
     */
    public static SVIntervalLocator of(final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary);
        final int numberOfContigs = dictionary.size();
        final int[] contigLengths = new int[numberOfContigs];
        final String[] contigNames = new String[numberOfContigs];
        for (int i = 0; i < numberOfContigs; i++) {
            final SAMSequenceRecord record = dictionary.getSequence(i);
            contigNames[i] = record.getSequenceName();
            contigLengths[i] = record.getSequenceLength();
        }

        return new SVIntervalLocator(contigNames, contigLengths);
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
        final int locStart = loc.getStart();
        final int locEnd = loc.getEnd();
        final int contigLength = contigLengths[index];
        Utils.validate(contigLength >= locEnd, "the input locatable runs beyond the end of the contig");
        Utils.validate(locStart >= 1 && locStart <= locEnd, "the input locatable start is of bounds");
        final int start = Math.max(1, loc.getStart() - padding);
        final int end = Math.min(contigLength + 1, locEnd + 1 + padding);
        return new SVInterval(index, start, end);
    }

    public <L extends Locatable, V> Collector<L, ?, SVIntervalTree<V> > toSVIntervalTree(final Function<L, V> toValue) {
        return Collector.of(SVIntervalTree<V>::new,
                            (tree, loc) -> tree.put(toSVInterval(loc), toValue.apply(loc)),
                            (tree1, tree2) -> { tree1.putAll(tree2); return tree1; });
    }

    public <L extends Locatable> L toLocatable(final SVInterval svInterval, final LocatableFactory<L> factory) {
        return toLocatable(svInterval, 0, factory);
    }

    public <L extends Locatable> L toLocatable(final SVInterval svInterval, final int padding, final LocatableFactory<L> factory) {
        Utils.nonNull(factory);
        Utils.nonNull(svInterval);
        ParamUtils.isPositiveOrZero(padding, "padding must be 0 or greater");
        final int index = svInterval.getContig();
        Utils.validateArg(index >= 0 && index < contigByIndex.length, "the input sv-interval index points to unexistent index");
        final int intervalStart = svInterval.getStart();
        final int intervalEnd = svInterval.getEnd() - 1;
        Utils.validate(intervalEnd <= contigLengths[index], "the input sv-interval goes beyond the contig end");
        return factory.create(contigByIndex[svInterval.getContig()],
                Math.max(1, intervalStart - padding), Math.min(intervalEnd + padding, contigLengths[index]));
    }

    public SimpleInterval toSimpleInterval(final SVInterval svInterval) {
        return toSimpleInterval(svInterval, 0);
    }

    public SimpleInterval toSimpleInterval(final SVInterval svInterval, final int padding) {
        return toLocatable(svInterval, padding, SimpleInterval::new);
    }

    public static class Serializer extends com.esotericsoftware.kryo.Serializer<SVIntervalLocator> {

        @Override
        public void write(final Kryo kryo, final Output output, final SVIntervalLocator object) {
            final int numberOfContigs = object.contigByIndex.length;
            output.writeInt(numberOfContigs);
            for (int i = 0; i < numberOfContigs; i++) {
                output.writeString(object.contigByIndex[i]);
                output.writeInt(object.contigLengths[i]);
            }
        }

        @Override
        public SVIntervalLocator read(Kryo kryo, Input input, Class<SVIntervalLocator> type) {
            final int numberOfContigs = input.readInt();
            final String[] contigsByIndex = new String[numberOfContigs];
            final int[] contigLenghts = new int[numberOfContigs];
            for (int i = 0; i < numberOfContigs; i++) {
                contigsByIndex[i] = input.readString();
                contigLenghts[i] = input.readInt();
            }
            return new SVIntervalLocator(contigsByIndex, contigLenghts);
        }
    }

    public SAMSequenceDictionary asDictionary() {
        final SAMSequenceDictionary result = new SAMSequenceDictionary();
        for (int i = 0; i < contigLengths.length; i++) {
            result.addSequence(new SAMSequenceRecord(contigByIndex[i], contigLengths[i]));
        }
        return result;
    }
}
