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
    private final int repeats;
    private final byte[] unit;

    private DragstrLocus(final int chrIdx, final long start, final byte[] unit, final int repeats) {
        chromosomeIndex = chrIdx;
        this.start = start;
        this.unit = unit;
        this.repeats = repeats;
    }

    public static DragstrLocus make(final int chrIdx, final long start, final byte[] unit, final int repeatCount) {
        ParamUtils.isPositiveOrZero(chrIdx, "chromosome index");
        ParamUtils.isPositive(start, "start position");
        final byte[] unitCloned = Utils.nonNull(unit).clone();
        for (final byte base : unit) {
            if (!Nucleotide.decode(base).isStandard()) {
                throw new IllegalArgumentException("bad bases in " + new String(unit));
            }
        }
        ParamUtils.isPositive(repeatCount, "repeat count");
        return new DragstrLocus(chrIdx, start, unit.clone(), repeatCount);
    }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return start + unit.length * repeats - 1;
    }

    public int getPeriod() {
        return unit.length;
    }

    byte[] getUnitUnsafe() {
        return unit;
    }

    public byte[] getUnit() {
        return unit.clone();
    }

    public int getRepeats() {
        return repeats;
    }

    public static BinaryTableWriter<DragstrLocus> binaryWriter(final OutputStream out, final String path) {
        return new BinaryTableWriter<DragstrLocus>(out, path) {
            @Override
            protected void writeRecord(final DragstrLocus record, final DataOutput output) throws IOException {
                output.writeShort(record.chromosomeIndex);
                output.writeLong(record.start);
                output.writeByte(record.unit.length);
                output.write(record.unit);
                output.writeByte(record.repeats);
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
                final int unitLength =  input.readUnsignedByte();
                final byte[] unit = new byte[unitLength];
                input.readFully(unit);
                final int repeats = input.readUnsignedByte();
                return new DragstrLocus(chrIdx, start, unit, repeats);
            }
        };
    }

    public static TableWriter<DragstrLocus> textWriter(final OutputStream out) throws IOException {
        return new TableWriter<DragstrLocus>(new OutputStreamWriter(out),
                TableColumnCollection.make("chridx", "start", "end", "length", "repeats", "unit")) {

            @Override
            protected void composeLine(final DragstrLocus record, final DataLine dataLine) {
                dataLine.append(record.chromosomeIndex)
                        .append(record.start)
                        .append(record.start + record.unit.length * record.repeats - 1)
                        .append(record.unit.length)
                        .append(record.repeats)
                        .append(new String(record.unit));
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
                final long end = dataLine.getLong("end");
                final int length = dataLine.getInt("length");
                final int repeats = dataLine.getInt("repeats");
                final byte[] unit = dataLine.get("unit").getBytes();
                final DragstrLocus site = new DragstrLocus(chridx, start, unit, repeats);
                if (site.getEnd() != end) {
                    throw formatException("end is not what is expected");
                } else if (site.getPeriod() != length) {
                    throw formatException("unit length is not what is expected");
                }
                return site;
            }
        };
    }

    public static TableReader<DragstrLocus> textReader(final InputStream in) throws IOException {
        return new TableReader<DragstrLocus>(new InputStreamReader(in)) {

            @Override
            protected DragstrLocus createRecord(final DataLine dataLine) {

                final int chridx = dataLine.getInt("chridx");
                final long start = dataLine.getLong("start");
                final long end = dataLine.getLong("end");
                final int length = dataLine.getInt("length");
                final int repeats = dataLine.getInt("repeats");
                final byte[] unit = dataLine.get("unit").getBytes();
                final DragstrLocus site = new DragstrLocus(chridx, start, unit, repeats);
                if (site.getEnd() != end) {
                    throw formatException("end is not what is expected");
                } else if (site.getPeriod() != length) {
                    throw formatException("unit length is not what is expected");
                }
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
        return ((((((chromosomeIndex * 31) + repeats) * 31) + (int) start) * 31 + Arrays.hashCode(unit)));
    }

    public boolean equals(final DragstrLocus other) {
        return other == this || (other.chromosomeIndex == this.chromosomeIndex &&
                other.repeats == this.repeats && other.start == this.start && equalUnits(this.unit, other.unit));
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
