package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.SimpleFormatter;

/**
 * Created by valentin on 10/21/17.
 */
@DefaultSerializer(AbstractSVHaplotype.Serializer.class)
public abstract class AbstractSVHaplotype implements SVHaplotype {

    protected final String name;
    protected final List<AlignmentInterval> intervals;
    protected final String variantId;
    protected final SimpleInterval variantLocation;

    @Override
    public List<AlignmentInterval> getReferenceAlignmentIntervals() {
        return intervals;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public String getVariantId() { return variantId; }

    @Override
    public SimpleInterval getVariantLocation() { return variantLocation; }

    protected AbstractSVHaplotype(final String name, final List<AlignmentInterval> intervals, final String variantId, final SimpleInterval variantLocation) {
        this.variantId = Utils.nonNull(variantId);
        this.variantLocation = Utils.nonNull(variantLocation).getStartInterval();
        this.name = Utils.nonNull(name);
        Utils.nonNull(intervals, "the input intervals cannot be null");
        ParamUtils.doesNotContainNulls(intervals, "the input intervals cannot contain nulls");
        //Utils.validateArg(!intervals.isEmpty(), "the input intervals cannot be empty");
        this.intervals = intervals.isEmpty() ? Collections.emptyList() : (intervals.size() == 1
                ? Collections.singletonList(intervals.get(0)) : Collections.unmodifiableList(new ArrayList<>(intervals)));
    }

    protected AbstractSVHaplotype(final Kryo kryo, final Input input) {
        this.name = input.readString();
        this.variantId = input.readString();
        final String variantContig = input.readString();
        final int variantStart = input.readInt();
        this.variantLocation = new SimpleInterval(variantContig, variantStart, variantStart);
        final int numberOfIntervals = input.readInt();
        final List<AlignmentInterval> intervals = new ArrayList<>(numberOfIntervals);
        for (int i = 0; i < numberOfIntervals; i++) {
            intervals.add(kryo.readObject(input, AlignmentInterval.class));
        }
        this.intervals = intervals.size() == 1 ? Collections.singletonList(intervals.get(0)) : Collections.unmodifiableList(new ArrayList<>(intervals));
    }

    /**
     * Copies haplotypes bases into an array assuming that the argument values passed are correct.
     *
     * @param offset
     * @param dest
     * @param destOffset
     * @param length
     */
    protected abstract void unsafeCopyBases(final int offset, final byte[] dest, final int destOffset, final int length);

    /**
     * Copies bases from the haplotype into an array.
     * @param offset
     * @param whereTo
     * @param destOffset
     * @param length
     * @throws IllegalArgumentException if {@code whereTo} is {@code null} or the indeces and length passed are not correct.
     */
    public void copyBases(final int offset, final byte[] whereTo, final int destOffset, final int length) {
        ParamUtils.isPositiveOrZero(length, "length must be positive or zero");
        ParamUtils.isPositiveOrZero(offset, "dest offset");
        ParamUtils.isPositiveOrZero(destOffset, "the offset must be positive or zero");
        Utils.validateArg(destOffset + length <= whereTo.length, "the to index must be less than the length of the haplotype");
        Utils.validateArg(offset + length <= getLength(), "the from index cannot be larger than the to index");
        Utils.nonNull(whereTo);
        unsafeCopyBases(offset, whereTo, destOffset, length);
    }

    @Override
    public int insertSize(final int start, final int end) {
        return end - start + 1;
    }

    public static class Serializer<S extends AbstractSVHaplotype> extends com.esotericsoftware.kryo.Serializer<S> {

        @Override
        public void write(Kryo kryo, Output output, S object) {
            output.writeString(object.name);
            output.writeString(object.variantId);
            output.writeString(object.variantLocation.getContig());
            output.writeInt(object.variantLocation.getStart());
            output.writeInt(object.intervals.size());
            for (final AlignmentInterval interval : object.intervals) {
                kryo.writeObject(output, interval);
            }
        }

        @Override
        public S read(Kryo kryo, Input input, Class<S> type) {
            throw new UnsupportedOperationException("read operation must be called at the specific class serializer");
        }
    }
}
