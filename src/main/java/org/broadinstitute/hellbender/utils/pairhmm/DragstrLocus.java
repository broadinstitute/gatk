package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.*;
import java.util.Arrays;

public class DragstrLocus {

    private final int chromosomeIndex;
    private final long start;
    private final byte period;
    private final short length;
    private final int mask;

    private DragstrLocus(final int chrIdx, final long start, final byte period, final short length, final int mask) {
        chromosomeIndex = chrIdx;
        this.start = start;
        this.period = period;
        this.length = length;
        this.mask = mask;
    }

    public static DragstrLocus make(final int chrIdx, final long start, final byte period, final short length, final int mask) {
        ParamUtils.isPositiveOrZero(chrIdx, "chromosome index");
        ParamUtils.isPositive(start, "start position");
        ParamUtils.isPositive(period, "period");
        ParamUtils.isPositive(length, "length");
        return new DragstrLocus(chrIdx, start, period, length, mask);
    }

    public long getMask() { return mask; }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return start + length - 1;
    }

    public int getPeriod() {
        return period;
    }

    public int getRepeats() {
        return length / period;
    }

    public static BinaryTableWriter<DragstrLocus> binaryWriter(final OutputStream out, final String path) {
        return new BinaryTableWriter<DragstrLocus>(out, path) {
            @Override
            protected void writeRecord(final DragstrLocus record, final DataOutput output) throws IOException {
                output.writeShort(record.chromosomeIndex);
                output.writeLong(record.start);
                output.writeByte(record.period);
                output.writeShort(record.length);
                output.writeInt(record.mask);
            }
        };
    }

    public static BinaryTableWriter<DragstrLocus> binaryWriter(final File out)
        throws FileNotFoundException
    {
        return binaryWriter(new FileOutputStream(out), out.toString());
    }

    public static BinaryTableReader<DragstrLocus> binaryReader(final InputStream in) {
        return new BinaryTableReader<DragstrLocus>(new DataInputStream(in)) {
            @Override
            protected DragstrLocus readRecord(final PushbackDataInput input) throws IOException {
                final int chrIdx = input.readUnsignedShort();
                final long start = input.readLong();
                final byte period = input.readByte();
                final short length = input.readShort();
                final int mask = input.readInt();
                return new DragstrLocus(chrIdx, start, period, length, mask);
            }
        };
    }

    public static TableWriter<DragstrLocus> textWriter(final OutputStream out) throws IOException {
        return new TableWriter<DragstrLocus>(new OutputStreamWriter(out),
                TableColumnCollection.make("chridx", "start", "period", "length", "mask")) {

            @Override
            protected void composeLine(final DragstrLocus record, final DataLine dataLine) {
                dataLine.append(record.chromosomeIndex)
                        .append(record.start)
                        .append(record.period)
                        .append(record.length)
                        .append(record.mask);
            }
        };
    }

    public static TableReader<DragstrLocus> textReader(final InputStream in, final SAMSequenceDictionary dictionary) throws IOException {
        return new TableReader<DragstrLocus>(new InputStreamReader(in)) {

            @Override
            protected DragstrLocus createRecord(final DataLine dataLine) {

                final String chr = dataLine.get("chr");
                final SAMSequenceRecord seq = dictionary.getSequence(chr);
                final int chridx = seq.getSequenceIndex();
                final long start = dataLine.getLong("start");
                final byte period = dataLine.getByte("period");
                final short length = (short) dataLine.getInt("length");
                final int mask = dataLine.getInt("mask");
                final DragstrLocus site = new DragstrLocus(chridx, start, period, length, mask);
                return site;
            }
        };
    }



    @Override
    public boolean equals(final Object other) {
        return (other instanceof DragstrLocus) && equals((DragstrLocus)other);
    }

    @Override
    public int hashCode() {
        return ((((((chromosomeIndex * 31) + length) * 31) + (int) start) * 31 + period));
    }

    public boolean equals(final DragstrLocus other) {
        return other == this || (other.chromosomeIndex == this.chromosomeIndex &&
                other.length == length && other.start == this.start && this.period == other.period);
    }

    private static boolean equalUnits(byte[] a, byte[] b) {
        if (a == b) {
            return true;
        } else if (a.length != b.length) {
            return false;
        } else {
            for (int i = 0; i < a.length; i++) {
                if (!Nucleotide.same(a[i],b[i])) {
                    return false;
                }
            }
            return true;
        }
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

    public SimpleInterval  getRepeatInterval(final SAMSequenceDictionary dictionary, final int leftMargin, final int rightMargin) {
        if (leftMargin == 0 && rightMargin == 0) {
            return new SimpleInterval(dictionary.getSequence(chromosomeIndex).getSequenceName(), (int) start, (int) getEnd());
        } else {
            final SAMSequenceRecord record = dictionary.getSequence(chromosomeIndex);
            final int start = Math.max(1, (int) this.start - leftMargin);
            final int end = Math.min(record.getSequenceLength(), (int) getEnd() + rightMargin);
            return new SimpleInterval(record.getSequenceName(), start, end);
        }
    }
}
