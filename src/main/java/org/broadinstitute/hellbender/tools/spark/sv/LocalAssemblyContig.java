package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.esotericsoftware.kryo.serializers.CollectionSerializer;
import com.esotericsoftware.kryo.serializers.DefaultArraySerializers;
import com.esotericsoftware.kryo.serializers.FieldSerializer;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Represents a collection of relevant information on one locally assembled contig from short reads.
 */
@DefaultSerializer(LocalAssemblyContig.Serializer.class)
class LocalAssemblyContig {

    /**
     * This contig was assembled based on FASTQ records collected by {@link FindBadGenomicKmersSpark} and numbered by it.
     */
    final long assemblyID;

    /**
     * This together with {@link #assemblyID} completely determines which contig we are talking about.
     */
    final String contigID;

    final String seq;

    /**
     * Read count for the contigs.
     * May be {@code null} if assembler doesn't provide this feature (e.g. SGA).
     */
    final int readCountSupport;
    static final int NA_READ_COUNT_SUPPORT = -1;

    /**
     * Per base coverage for the assembled contig
     * May be {@code null} if assembler doesn't provide this feature (e.g. SGA).
     */
    @FieldSerializer.Bind(DefaultArraySerializers.ByteArraySerializer.class)
    final byte[] perBaseCoverage;

    @CollectionSerializer.BindCollection(elementClass = AlignmentRegion.class, elementSerializer = AlignmentRegion.Serializer.class, elementsCanBeNull = true)
    ArrayList<AlignmentRegion> alignmentRecords; // using ArrayList for serialization reason

    /**
     * Construct a contig with assembly id, contig/vertex id, bases, but without number of reads generating the contig, or qual.
     */
    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq){
        this(assemblyID, contigID, seq, null, NA_READ_COUNT_SUPPORT, null);
    }

    //using ArrayList for serialization reason
    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq,
                        final ArrayList<AlignmentRegion> alignmentRegions){
        this(assemblyID, contigID, seq, alignmentRegions, NA_READ_COUNT_SUPPORT, null);
    }

    /**
     * Construct a contig with assembly id, contig/vertex id, bases, number of reads generating the contig, and qual.
     * Currently fermi-lite supports the last two features, not SGA.
     * Using ArrayList for serialization reason
     */
    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq,
                        final ArrayList<AlignmentRegion> alignmentRegions,
                        final int readCountSupport, final byte[] perBaseCoverage){
        this.assemblyID = assemblyID;
        this.contigID   = contigID;
        this.seq        = seq;
        this.readCountSupport = readCountSupport;

        this.perBaseCoverage = (perBaseCoverage == null) ? null : perBaseCoverage;
        this.alignmentRecords = (alignmentRegions == null) ? null : alignmentRegions;
    }

    void setAlignmentRegions(final Collection<AlignmentRegion> ars) {
        alignmentRecords = ars==null ? null : new ArrayList<>(ars);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        LocalAssemblyContig that = (LocalAssemblyContig) o;

        if (assemblyID != that.assemblyID) return false;
        if (!contigID.equals(that.contigID)) return false;
        return seq.equals(that.seq);
    }

    @Override
    public int hashCode() {
        int result = (int) (assemblyID ^ (assemblyID >>> 32));
        result = 31 * result + contigID.hashCode();
        result = 31 * result + seq.hashCode();
        return result;
    }

    ////////////////////////////////////////// Serialization stuff
    @SuppressWarnings("unchecked")
    protected LocalAssemblyContig(final Kryo kryo, final Input input){
        this.assemblyID = input.readLong();
        this.contigID = input.readString();
        this.seq = input.readString();
        this.readCountSupport = input.readInt();

        this.perBaseCoverage = kryo.readObjectOrNull(input, byte[].class); // TODO: how do we know the length?
        this.alignmentRecords = kryo.readObjectOrNull(input, ArrayList.class);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<LocalAssemblyContig>{
        @Override
        public LocalAssemblyContig read(final Kryo kryo, final Input input, final Class<LocalAssemblyContig> klass){
            return new LocalAssemblyContig(kryo, input);
        }

        @Override
        public void write(final Kryo kryo, final Output output, final LocalAssemblyContig contig){
            contig.serialize(kryo, output);
        }
    }

    protected void serialize(final Kryo kryo, final Output output){
        output.writeLong(assemblyID);
        output.writeString(contigID);
        output.writeString(seq);
        output.writeInt(readCountSupport);

        kryo.writeObjectOrNull(output, perBaseCoverage, byte[].class);
        kryo.writeObjectOrNull(output, alignmentRecords, ArrayList.class);
    }
}
