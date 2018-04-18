package org.broadinstitute.hellbender.utils.read.markduplicates;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

public abstract class MarkDuplicatesSparkData {
    protected final int partitionIndex;
    protected final String name;

    MarkDuplicatesSparkData(int partitionIndex, String name) {
        this.name = name;
        this.partitionIndex = partitionIndex;
    }

    public abstract Type getType();
    public abstract int getScore();
    protected abstract void serialize(Kryo kryo, Output output);
    public abstract int getUnclippedStartPosition();
    public abstract int getFirstRefIndex();
    public abstract int key(final SAMFileHeader header);


    public int getPartitionIndex(){
      return partitionIndex;
    }
    public String getName() {
      return name;
    }

    public enum Type {
        FRAGMENT, PAIR, PASSTHROUGH;
    }

    public static Passthrough getPassthrough(GATKRead read, int partitionIndex) {
        return new Passthrough(read, partitionIndex);
    }

    /**
     * Dummy class used for preserving reads that need to be marked as non-duplicate despite not wanting to perform any
     * processing on the reads. (eg. unmapped reads we don't want to process but must be non-duplicate marked)
     */
    public static final class Passthrough extends MarkDuplicatesSparkData {
        private final transient GATKRead read;

        Passthrough(GATKRead read, int partitionIndex) {
            super(partitionIndex, read.getName());

            this.read = read;
        }

        @Override
        public Type getType() {
            return Type.PASSTHROUGH;
        }
        // Dummy values for an empty class
        @Override
        public int getScore() { return 0; }
        @Override
        public int getUnclippedStartPosition() { return 0; }
        @Override
        public int getFirstRefIndex() { return 0; }
        @Override
        public int key(final SAMFileHeader header) {
            return ReadsKey.hashKeyForRead(read);
        }

        // Serialization support
        @Override
        protected void serialize(Kryo kryo, Output output) {
            output.writeInt(partitionIndex, true);
            output.writeAscii(name);
        }
        Passthrough(Kryo kryo, Input input) {
            super(input.readInt(true), input.readString());
            this.read = null;
        }
    }

    public static final class MarkDuplicatesSparkDataPassthroughSerializer extends com.esotericsoftware.kryo.Serializer<MarkDuplicatesSparkData.Passthrough> {
        @Override
        public void write(final Kryo kryo, final Output output, final MarkDuplicatesSparkData.Passthrough pair ) {
            pair.serialize(kryo, output);
        }

        @Override
        public MarkDuplicatesSparkData.Passthrough read(final Kryo kryo, final Input input, final Class<MarkDuplicatesSparkData.Passthrough> klass ) {
            return new MarkDuplicatesSparkData.Passthrough(kryo, input);
        }
    }
}
