package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;

import java.util.stream.LongStream;

/**
 * Contains tuple of counts associated with a target.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReadCountRecord implements Feature {

    private final Target target;

    // counts, one element per read count group.
    private final double[] counts;

    /**
     * Creates a new read-count record given the target and the counts.
     * <p>
     *     The resulting record has its own copy of the counts, therefore future changes in the input count
     *     array will not affect its state.
     * </p>
     * @param target the target.
     * @param counts the counts for the new record.
     * @throws IllegalArgumentException if either {@code target} or {@code counts} is {@code null}.
     */
    public ReadCountRecord(final Target target, final double[] counts) {
        this.target = Utils.nonNull(target);
        this.counts = Utils.nonNull(counts).clone();
    }

    /**
     * Creates a new read-count record given the target and the counts.
     * <p>
     *     The resulting record has its own copy of the counts, therefore future changes in the input count
     *     array will not affect its state.
     * </p>
     * @param target the target.
     * @param counts the counts for the new record.
     * @throws IllegalArgumentException if either {@code target} or {@code counts} is {@code null}.
     */
    public ReadCountRecord(final Target target, final long[] counts) {
        this.target = Utils.nonNull(target);
        this.counts = LongStream.of(Utils.nonNull(counts)).mapToDouble(l -> l).toArray();
    }

    /**
     * Queries this record target.
     *
     * @return never {@code null}.
     */
    public Target getTarget() {
        return target;
    }

    /**
     * Appends this instance read counts to a data-line object.
     *
     * @param dataLine the destination data-line.
     * @throws IllegalArgumentException if {@code dataLine is {@code null}}.
     * @throws IllegalStateException    if there is no enough room in {@code dataLine} from is current appending position
     *                                  to add al the record counts.
     */
    public void appendCountsTo(final DataLine dataLine) {
        Utils.nonNull(dataLine);
        dataLine.append(counts);
    }

    /**
     * copy this instance read counts to a double array from a particular position.
     * @param destination the array to copy the values to. This call will modify its content.
     * @param offset the position in the array for the first count in this record.
     *
     * @throws IllegalArgumentException if {@code destination} is {@code null}, {@code offset} is
     *   negative or {@code destination} is not long enough to accommodate for all
     *   the values of this record starting at position {@code offset}.
     */
    public void copyCountsTo(final double[] destination, final int offset) {
        Utils.nonNull(destination);
        Utils.validateArg(offset >= 0, () -> "the offset must be 0 or greater: " + offset);
        final int to = counts.length + offset;
        Utils.validateArg(to <= destination.length, () -> String.format("the destination array is not long enough: %d > %d", to, destination.length));
        System.arraycopy(counts, 0, destination, offset, counts.length);
    }

    /**
     * Return a copy of the counts in the record.
     *
     * @return never {@code null}.
     */
    public double[] getDoubleCounts() {
        return counts.clone();
    }

    /**
     * Number of counts in this read count record.
     *
     * @return 0 or greater.
     */
    public int size() {
        return counts.length;
    }

    /**
     * Return a count as a double.
     *
     * @param index the count index.
     * @return any double value.
     * @throws IllegalArgumentException if {@code index} is outside the valid range [0 .. {@link #size()} - 1].
     */
    public double getDouble(final int index) {
        Utils.validIndex(index, counts.length);
        return counts[index];
    }

    public SingleSampleRecord asSingleSampleRecord() {
        return new SingleSampleRecord(this);
    }

    @Override
    @Deprecated
    public String getChr() {
        return target.getContig();
    }

    @Override
    public String getContig() {
        return target.getContig();
    }

    @Override
    public int getStart() {
        return target.getStart();
    }

    @Override
    public int getEnd() {
        return target.getEnd();
    }

    public static class SingleSampleRecord extends ReadCountRecord {
        public SingleSampleRecord(final Target target, final double[] counts) {
            super(target, counts);
            Utils.validateArg(counts.length == 1, "SingleSampleRecord must have exactly one count.");
        }

        public SingleSampleRecord(final Target target, final double count) {
            this(target, new double[] {count});
        }

        public SingleSampleRecord(final ReadCountRecord record) {
            this(record.target, record.counts);
        }

        public double getCount() {
            return getDouble(0);
        }
     }
}
