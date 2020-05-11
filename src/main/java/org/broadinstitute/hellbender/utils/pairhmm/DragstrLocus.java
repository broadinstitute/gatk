package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import it.unimi.dsi.fastutil.ints.*;
import org.apache.hadoop.io.IOUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
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

public class DragstrLocus {

    private final int chromosomeIndex;
    private final long start;
    private final byte period;
    private final short length;
    private final int mask;

    private static final int INDEX_BYTE_INTERVAL = 1 << 16; // every 64KB

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

    public int getChromosomeIndex() { return chromosomeIndex; }

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

    @FunctionalInterface
    public interface WriteAction {
        void write(final DragstrLocus locus, final DataOutput dest) throws IOException;
    }


    public static BinaryTableWriter<DragstrLocus> binaryWriter(final OutputStream out, final OutputStream indexOut, final String path) {
        return binaryWriter(out, indexOut, path, (record, output) -> {
            output.writeInt(record.mask);
            output.writeShort(record.chromosomeIndex);
            output.writeLong(record.start);
            output.writeByte(record.period);
            output.writeShort(record.length);
        });
    }

    public static BinaryTableWriter<DragstrLocus> dragenWriter(final OutputStream out, final OutputStream indexOut, final String path) {
        final ByteBuffer buffer = ByteBuffer.allocate(Integer.BYTES * 3 + Short.BYTES + 2 * Byte.BYTES);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        return binaryWriter(out, indexOut, path, (record, output) -> {
            buffer.clear();
            buffer.putInt(record.mask);
            buffer.putInt(record.chromosomeIndex);
            buffer.putInt((int) record.start - 1);  
            buffer.putShort(record.length);
            buffer.put(record.period);
            buffer.put((byte) Math.min(255, record.getRepeats()));
            output.write(buffer.array());
        });
    }


    private static BinaryTableWriter<DragstrLocus> binaryWriter(final OutputStream out, final OutputStream indexOut, final String path, final WriteAction wa) {
        return new BinaryTableWriter<DragstrLocus>(out, path) {

            private DataOutputStream indexDataOutputStream = indexOut != null ? new DataOutputStream(indexOut) : null;
            private Int2LongMap chromosomeOffsets = new Int2LongArrayMap();
            private int lastChromosomeIndex = -1;
            private long lastEntryOffset = 0;

            @Override
            protected void writeRecord(final DragstrLocus record, final DataOutput output) throws IOException {
                if (indexOut != null) {
                    outputIndexWhenApplies(record);
                }
                wa.write(record, output);
            }

            private void outputIndexWhenApplies(DragstrLocus record) throws IOException {
                final long offset = offset();
                if (lastChromosomeIndex != record.chromosomeIndex) {
                    if (chromosomeOffsets.containsKey(record.chromosomeIndex)) {
                        throw new IllegalStateException("cannot index the output when not sorted; chromosome idex " + record.chromosomeIndex + " apears in more than one piece ");
                    }
                    dataOut.flush();
                    chromosomeOffsets.put(lastChromosomeIndex = record.chromosomeIndex, offset);
                    indexDataOutputStream.writeInt(lastChromosomeIndex);
                    indexDataOutputStream.writeInt((int) record.start);
                    indexDataOutputStream.writeLong(offset);
                    lastEntryOffset = offset;
                } else if ((offset - lastEntryOffset) >= INDEX_BYTE_INTERVAL) {
                    dataOut.flush();
                    indexDataOutputStream.writeInt(lastChromosomeIndex);
                    indexDataOutputStream.writeInt((int) record.start);
                    indexDataOutputStream.writeLong(offset);
                    lastEntryOffset = offset;
                }
            }

            @Override
            public void close() throws IOException {
                super.close();
                if (indexOut != null) indexDataOutputStream.close();
            }
        };
    }

    public static BinaryTableWriter<DragstrLocus> binaryWriter(final File out) throws FileNotFoundException {
        return binaryWriter(out, null);
    }

    public static BinaryTableWriter<DragstrLocus> binaryWriter(final File out, final File indexFile)
        throws FileNotFoundException
    {
        return binaryWriter(new FileOutputStream(out), indexFile != null ? new FileOutputStream(indexFile) : new IOUtils.NullOutputStream(), out.toString());
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



    public static BinaryTableReader<DragstrLocus> binaryReader(final String path, final BinaryTableIndex index, final int chrIdx, final int start, final int end)
       throws IOException {

        final long offset = index.offset(chrIdx, start, end);
        if (offset < 0) {
            return BinaryTableReader.emptyReader();
        }
        final InputStream in = BucketUtils.openFile(path);
        in.skip(offset);

        return new BinaryTableReader<DragstrLocus>(new DataInputStream(in)) {

            @Override
            protected DragstrLocus readRecord(final PushbackDataInput input) throws IOException {
                while (!input.eof()) {
                    final int c = input.readUnsignedShort();
                    final long s = input.readLong();
                    final byte p = input.readByte();
                    final short l = input.readShort();
                    final int m = input.readInt();
                    final long e = s + l - 1;
                    if (chrIdx != c) {
                        return null;
                    } else if (s > end) {
                        return null;
                    } else if (s >= start) {
                        return new DragstrLocus(c, s, p, l, m);
                    }
                }
                return null;
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

    public static class BinaryTableIndex {

        private static class Entry {
            public final int chrIdx;
            public final int pos;
            public final long offset;

            private Entry(final int chrIdx, final int pos, final long offset) {
                this.chrIdx = chrIdx;
                this.pos = pos;
                this.offset = offset;
            }

            public static Entry of(final int chrIdx, final int pos, final long offset) {
                return new Entry(chrIdx, pos, offset);
            }

        }

        private Int2ObjectMap<List<Entry>> entriesByChrIdx;

        private BinaryTableIndex(final Int2ObjectMap<List<Entry>> entries) {
            entriesByChrIdx = entries;
        }

        public long offset(final int chrIdx, final int start, final int end) {
            final List<Entry> chrEntries = entriesByChrIdx.get(chrIdx);
            if (chrEntries == null || chrEntries.isEmpty()) {
                return -1;
            } else if (chrEntries.get(0).pos > end) {
                return -1;
            } else {
                int i = 0, j = chrEntries.size() - 1;
                while (i < j) {
                    int k = (i + j) / 2;
                    final Entry candidate = chrEntries.get(k);
                    if (candidate.pos < start) {
                        i = Math.min(k + 1, j);
                    } else if (candidate.pos > start) {
                        j = Math.max(k - 1, i);
                    } else {
                        return candidate.offset;
                    }
                }
                while (i > 0 && chrEntries.get(i).pos > start) {
                    i--;
                }
                return chrEntries.get(i).offset;
            }
        }

        public static BinaryTableIndex load(final String path) throws IOException {
            return load(BucketUtils.openFile(path));
        }

        public static BinaryTableIndex load(final InputStream inputStream) throws IOException {

            final BinaryTableReader<Entry> entryReader = new BinaryTableReader<Entry>(new DataInputStream(inputStream)) {
                @Override
                protected Entry readRecord(PushbackDataInput input) throws IOException {
                    final int chrIdx = input.readInt();
                    final int pos = input.readInt();
                    final long offset = input.readLong();
                    return Entry.of(chrIdx, pos, offset);
                }
            };
            final Int2ObjectMap<List<Entry>> entriesByChrIdx = new Int2ObjectArrayMap<>();
            entryReader.stream()
                       .forEach(entry -> {
                           List<Entry> entries = entriesByChrIdx.get(entry.chrIdx);
                           if (entries == null) {
                               entriesByChrIdx.put(entry.chrIdx, entries = new ArrayList<>());
                           }
                           entries.add(entry);
                       });

            return new BinaryTableIndex(entriesByChrIdx);
        }
    }
}
