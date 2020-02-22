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

    private static final int INDEX_BYTE_INTERVAL = 1 << 16; // every 64KB

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


    private static BinaryTableWriter<DragstrLocus> binaryWriter(final OutputStream out, final OutputStream indexOut, final String path) {
        return binaryWriter(out, indexOut, path, (record, output) -> {
            output.writeInt(record.chromosomeIndex);
            output.writeLong(record.start);
            output.writeByte(record.period);
            output.writeShort(record.length);
            output.writeLong(record.mask);
        });
    }

    /**
     * Generates the str_table.bin format produced by DRAGEN. Used for debugging purposes. Please keep around for now.
     */
    @SuppressWarnings("unused")
    public static BinaryTableWriter<DragstrLocus> dragenWriter(final OutputStream out, final OutputStream indexOut, final String path) {
        final ByteBuffer buffer = ByteBuffer.allocate(Integer.BYTES * 3 + Short.BYTES + 2 * Byte.BYTES);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        return binaryWriter(out, indexOut, path, (record, output) -> {
            buffer.clear();
            buffer.putInt((int) record.mask);
            buffer.putInt(record.chromosomeIndex);
            buffer.putInt((int) record.start - 1);
            buffer.putShort(record.length);
            buffer.put(record.period);
            buffer.put((byte) Math.min(20, record.getRepeats())); // DRAGEN caps repeat lengths to 20.
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

            private void outputIndexWhenApplies(final DragstrLocus record) throws IOException {
                final long offset = offset();
                if (lastChromosomeIndex != record.chromosomeIndex) {
                    if (chromosomeOffsets.containsKey(record.chromosomeIndex)) {
                        throw new IllegalStateException("cannot index the output when not sorted; chromosome index " + record.chromosomeIndex + " appears in more than one piece ");
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


    public static BinaryTableReader<DragstrLocus> binaryReader(final File file) throws FileNotFoundException {
        return binaryReader(new FileInputStream(file));
    }

    public static BinaryTableReader<DragstrLocus> binaryReader(final InputStream in) {
        return new BinaryTableReader<DragstrLocus>(in, null) {
            @Override
            protected DragstrLocus readRecord(final DataInput input) throws IOException {
                final int chrIdx = input.readInt();
                final long start = input.readLong();
                final byte period = input.readByte();
                final short length = input.readShort();
                final long mask = input.readLong();
                return DragstrLocus.make(chrIdx, start, period, length, mask);
            }
        };
    }

    /**
     * Returns loci whose start base is located within an particular interval.
     * @param path path to the file containing the dragstr-loci.
     * @param index the index for that file pre-loaded in memory.
     * @param chrIdx the target interval contig index.
     * @param start the first base of the target interval 1-based.
     * @param end the last base of the target interval (inclusive).
     * @return never {@code null} but perhaps a read that returns not records.
     * @throws IOException in case of an underlying IO issue.
     */
    public static BinaryTableReader<DragstrLocus> binaryReader(final String path, final BinaryTableIndex index,
                                                               final int chrIdx, final int start, final int end)
       throws IOException {

        final long offset = index.offset(chrIdx, start, end);
        if (offset < 0) {
            return BinaryTableReader.emptyReader();
        }
        final InputStream in = BucketUtils.openFile(path);
        if (in.skip(offset) != offset) {
            throw new IOException("failed to skip the requested number of bytes");
        }

        return new BinaryTableReader<DragstrLocus>(in, null) {

            @Override
            protected DragstrLocus readRecord(final DataInput input) throws IOException {

                while (true) {
                    final int c = input.readInt();
                    final long s = input.readLong();
                    final byte p = input.readByte();
                    final short l = input.readShort();
                    final long m = input.readLong();
                    if (chrIdx != c) { // always we have an entry in the index for each chromosome
                                       // so we should not any chridx that is not the intervals is end of the line.
                        return null;
                    } else if (s > end) {
                        return null;
                    } else if (s >= start) {
                        return new DragstrLocus(c, s, p, l, m);
                    }
                    // notice that eventually we hit the end of the stream or another chromosome.
                    // or we go beyond the requested interval, so this is not going to loop forever.
                }
            }
        };
    }

    /**
     * Returns a tab separated text format writer.
     * @param out the output stream where to write the data to.
     * @param dictionary the dictionary for the reference this loci refer to.
     * @return never {@code null}.
     * @throws IOException iff there is any low-level issue creating the writer.
     */
    public static TableWriter<DragstrLocus> textWriter(final OutputStream out, final SAMSequenceDictionary dictionary) throws IOException {


        return new TableWriter<DragstrLocus>(new OutputStreamWriter(out),
                TableColumnCollection.make("chridx", "chrid", "start", "end", "period", "mask", "mask_bin", "length_bp", "length_rp")) {

            @Override
            protected void composeLine(final DragstrLocus record, final DataLine dataLine) {
                dataLine.append(record.chromosomeIndex)
                        .append(dictionary.getSequence(record.chromosomeIndex).getSequenceName())
                        .append(record.start)
                        .append(record.start + record.length - 1)
                        .append(record.period)
                        .append(record.mask)
                        .append(Long.toBinaryString(record.mask))
                        .append(record.length)
                        .append(record.length / record.period);
            }
        };
    }

    /**
     * Reads in DragstrLocus instances from a stream that has the same content as the one generated using the {@link #textWriter}.
     */
    static TableReader<DragstrLocus> textReader(final InputStream in, final SAMSequenceDictionary dictionary) throws IOException {
        return new TableReader<DragstrLocus>(new InputStreamReader(in)) {

            @Override
            protected DragstrLocus createRecord(final DataLine dataLine) {

                final String chr = dataLine.get("chrid");
                final SAMSequenceRecord seq = dictionary.getSequence(chr);
                final int chridx = seq.getSequenceIndex();
                final long start = dataLine.getLong("start");
                final byte period = dataLine.getByte("period");
                final short length = (short) dataLine.getInt("length");
                final int mask = dataLine.getInt("mask");
                return new DragstrLocus(chridx, start, period, length, mask);
            }
        };
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

        public static BinaryTableIndex load(final String path) {
            return load(BucketUtils.openFile(path));
        }

        public static BinaryTableIndex load(final InputStream inputStream) {

            final BinaryTableReader<Entry> entryReader = new BinaryTableReader<Entry>(inputStream, null) {
                @Override
                protected Entry readRecord(final DataInput input) throws IOException {
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
