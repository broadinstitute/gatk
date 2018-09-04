package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;

import java.util.ArrayList;
import java.util.List;

/**
 * Holding necessary information about a local assembly for use in SV discovery.
 */
@DefaultSerializer(AlignedAssembly.Serializer.class)
public final class AlignedAssembly {

    public final int assemblyId;

    public final List<AlignedContig> alignedContigs;

    public AlignedAssembly(final int assemblyId, final List<AlignedContig> alignedContigs) {
        this.assemblyId = assemblyId;
        this.alignedContigs = alignedContigs;
    }

    @VisibleForTesting
    private AlignedAssembly(final Kryo kryo, final Input input) {
        this.assemblyId = input.readInt();

        final int nContigs = input.readInt();
        alignedContigs = new ArrayList<>(nContigs);
        for(int contigIdx = 0; contigIdx < nContigs; ++contigIdx) {
            alignedContigs.add(new AlignedContig(kryo, input));
        }
    }

    @VisibleForTesting
    private void serialize(final Kryo kryo, final Output output) {
        output.writeInt(assemblyId);

        output.writeInt(alignedContigs.size());
        alignedContigs.forEach(alignedContig -> alignedContig.serialize(kryo, output));
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignedAssembly> {
        @Override
        public void write( final Kryo kryo, final Output output, final AlignedAssembly alignedAssembly){
            alignedAssembly.serialize(kryo, output);
        }

        @Override
        public AlignedAssembly read(final Kryo kryo, final Input input, final Class<AlignedAssembly> clazz ) {
            return new AlignedAssembly(kryo, input);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AlignedAssembly that = (AlignedAssembly) o;

        if (assemblyId != that.assemblyId) return false;
        return alignedContigs.equals(that.alignedContigs);
    }

    @Override
    public int hashCode() {
        int result = assemblyId;
        result = 31 * result + alignedContigs.hashCode();
        return result;
    }

}
