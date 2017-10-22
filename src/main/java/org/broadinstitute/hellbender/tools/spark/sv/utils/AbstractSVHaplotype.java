package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import org.apache.ivy.util.CollectionUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ops.impl.transforms.Abs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by valentin on 10/21/17.
 */
public abstract class AbstractSVHaplotype implements SVHaplotype {

    protected final String name;
    protected final List<AlignmentInterval> intervals;

    @Override
    public List<AlignmentInterval> getReferenceAlignmentIntervals() {
        return intervals;
    }

    @Override
    public String getName() {
        return name;
    }

    protected AbstractSVHaplotype(final String name, final List<AlignmentInterval> intervals) {
        this.name = Utils.nonNull(name);
        Utils.nonNull(intervals, "the input intervals cannot be null");
        ParamUtils.doesNotContainNulls(intervals, "the input intervals cannot contain nulls");
        Utils.validateArg(!intervals.isEmpty(), "the input intervals cannot be empty");
        this.intervals = intervals.size() == 1 ? Collections.singletonList(intervals.get(0)) : Collections.unmodifiableList(new ArrayList<>(intervals));
    }

    protected AbstractSVHaplotype(final Kryo kryo, final Input input) {
        this.name = input.readString();
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
}
