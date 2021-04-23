package org.broadinstitute.hellbender.tools.dragstr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import it.unimi.dsi.fastutil.ints.*;
import org.apache.hadoop.io.IOUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.BinaryTableReader;
import org.broadinstitute.hellbender.utils.BinaryTableWriter;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;

/**
 * Holds information about a locus on the reference that might be used to estimate the DRAGstr model parameters.
 */
public class DragstrLocus implements Comparable<DragstrLocus> {

    private final int chromosomeIndex;
    private final long start;
    private final byte period;
    private final short length;
    private final long mask;

    private DragstrLocus(final int chrIdx, final long start, final byte period, final short length, final long mask) {
        chromosomeIndex = chrIdx;
        this.start = start;
        this.period = period;
        this.length = length;
        this.mask = mask;
    }

    public static DragstrLocus make(final int chrIdx, final long start, final byte period, final short length, final long mask) {
        ParamUtils.isPositiveOrZero(chrIdx, "chromosome index");
        ParamUtils.isPositive(start, "start position");
        ParamUtils.isPositive(period, "period");
        ParamUtils.isPositive(length, "length");
        return new DragstrLocus(chrIdx, start, period, length, mask);
    }

    public int getChromosomeIndex() { return chromosomeIndex; }

    public long getMask() { return mask; }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return start + length - 1;
    }

    public short getLength() { return length; }

    public int getPeriod() {
        return period;
    }

    public int getRepeats() {
        return period == 0 ? 0 : length / period;
    }

    @Override
    public int compareTo(final DragstrLocus other) {
        if (other == this) {
            return 0;
        } else if (other == null) {
            return 1;
        } else {
            int cmp;
            if ((cmp = Integer.compare(this.chromosomeIndex, other.chromosomeIndex)) != 0) {
                return cmp;
            } else if ((cmp = Long.compare(this.start, other.start)) != 0) {
                return cmp;
            } else if ((cmp = Integer.compare(this.length, other.length)) != 0) {
                return cmp;
            } else {
                return Integer.compare(this.period, other.period);
            }
        }
    }

    @Override
    public boolean equals(final Object other) {
        return other instanceof DragstrLocus && compareTo((DragstrLocus) other) == 0;
    }

    @Override
    public int hashCode() {
        // in practice chridx and start position are enough as there should be rare to have collisions with only those two.
        return ((chromosomeIndex * 31) + (int) start) * 47;
    }

    @FunctionalInterface
    public interface WriteAction {
        void write(final DragstrLocus locus, final DataOutput dest) throws IOException;
    }


    public SimpleInterval getStartInterval(final SAMSequenceDictionary dictionary, final int margin) {
        if (margin == 0) {
            return new SimpleInterval(dictionary.getSequence(chromosomeIndex).getSequenceName(), (int) start, (int) start);
        } else {
            final SAMSequenceRecord record = dictionary.getSequence(chromosomeIndex);
            final int start = Math.max(1, (int) this.start - margin);
            final int end = Math.min(record.getSequenceLength(), (int) this.start + margin);
            return new SimpleInterval(record.getSequenceName(), start, end);
         }
    }

}
