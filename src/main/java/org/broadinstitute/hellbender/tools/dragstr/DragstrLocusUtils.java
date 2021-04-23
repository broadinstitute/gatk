package org.broadinstitute.hellbender.tools.dragstr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import it.unimi.dsi.fastutil.ints.Int2LongArrayMap;
import it.unimi.dsi.fastutil.ints.Int2LongMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import org.apache.hadoop.io.IOUtils;
import org.broadinstitute.hellbender.utils.BinaryTableReader;
import org.broadinstitute.hellbender.utils.BinaryTableWriter;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;

public class DragstrLocusUtils {
    private static final int INDEX_BYTE_INTERVAL = 1 << 16; // every 64KB

    private static BinaryTableWriter<DragstrLocus> binaryWriter(final OutputStream out, final OutputStream indexOut, final String path) {
        return binaryWriter(out, indexOut, path, (record, output) -> {
            output.writeInt(record.getChromosomeIndex());
            output.writeLong(record.getStart());
            output.writeByte(record.getPeriod());
            output.writeShort(record.getLength());
            output.writeLong(record.getMask());
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
            buffer.putInt((int) record.getMask());
            buffer.putInt(record.getChromosomeIndex());
            buffer.putInt((int) record.getStart() - 1);
            buffer.putShort(record.getLength());
            buffer.put((byte) record.getPeriod());
            buffer.put((byte) Math.min(DragstrHyperParameters.DEFAULT_MAX_PERIOD, record.getRepeats())); // DRAGEN caps repeat lengths to 20, the default max period.
            output.write(buffer.array());
        });
    }

    private static BinaryTableWriter<DragstrLocus> binaryWriter(final OutputStream out, final OutputStream indexOut, final String path, final DragstrLocus.WriteAction wa) {
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
                if (lastChromosomeIndex != record.getChromosomeIndex()) {
                    if (chromosomeOffsets.containsKey(record.getChromosomeIndex())) {
                        throw new IllegalStateException("cannot index the output when not sorted; chromosome index " + record.getChromosomeIndex() + " appears in more than one piece ");
                    }
                    dataOut.flush();
                    chromosomeOffsets.put(lastChromosomeIndex = record.getChromosomeIndex(), offset);
                    indexDataOutputStream.writeInt(lastChromosomeIndex);
                    indexDataOutputStream.writeInt((int) record.getStart());
                    indexDataOutputStream.writeLong(offset);
                    lastEntryOffset = offset;
                } else if ((offset - lastEntryOffset) >= INDEX_BYTE_INTERVAL) {
                    dataOut.flush();
                    indexDataOutputStream.writeInt(lastChromosomeIndex);
                    indexDataOutputStream.writeInt((int) record.getStart());
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
                        return DragstrLocus.make(c, s, p, l, m);
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
                dataLine.append(record.getChromosomeIndex())
                        .append(dictionary.getSequence(record.getChromosomeIndex()).getSequenceName())
                        .append(record.getStart())
                        .append(record.getStart() + record.getLength() - 1)
                        .append(record.getPeriod())
                        .append(record.getMask())
                        .append(Long.toBinaryString(record.getMask()))
                        .append(record.getLength())
                        .append(record.getLength() / record.getPeriod());
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
                return DragstrLocus.make(chridx, start, period, length, mask);
            }
        };
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
