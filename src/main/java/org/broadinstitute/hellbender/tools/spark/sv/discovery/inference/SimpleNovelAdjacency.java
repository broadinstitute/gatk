package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

// TODO: 1/19/18 step 1 towards the comprehensive view for Sv type inference
@DefaultSerializer(SimpleNovelAdjacency.Serializer.class)
public final class SimpleNovelAdjacency {

    private static final NovelAdjacencyAndInferredAltHaptype.Serializer narlSerializer = new NovelAdjacencyAndInferredAltHaptype.Serializer();
    private static final ChimericAlignment.Serializer caSerializer = new ChimericAlignment.Serializer();

    private final NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype;
    private final List<ChimericAlignment> alignmentEvidence;

    public NovelAdjacencyReferenceLocations getNovelAdjacencyReferenceLocations() {
        return novelAdjacencyAndInferredAltHaptype.novelAdjacencyReferenceLocations;
    }
    public byte[] getAltHaplotypeSequence() {
        return novelAdjacencyAndInferredAltHaptype.altHaplotypeSequence;
    }
    public List<ChimericAlignment> getAlignmentEvidence() {
        return alignmentEvidence;
    }

    public SimpleNovelAdjacency(final NovelAdjacencyAndInferredAltHaptype novelAdjacencyReferenceLocations,
                                final List<ChimericAlignment> alignmentEvidence) {
        this.novelAdjacencyAndInferredAltHaptype = Utils.nonNull( novelAdjacencyReferenceLocations );
        this.alignmentEvidence = Utils.nonNull( alignmentEvidence );
    }

    private SimpleNovelAdjacency(final Kryo kryo, final Input input) {
        novelAdjacencyAndInferredAltHaptype = narlSerializer.read(kryo, input, NovelAdjacencyAndInferredAltHaptype.class);
        final int evidenceCount = input.readInt();
        alignmentEvidence = new ArrayList<>(evidenceCount);
        for (int i = 0; i < evidenceCount; ++i) {
            alignmentEvidence.add(caSerializer.read(kryo, input, ChimericAlignment.class));
        }
    }

    private void serialize(final Kryo kryo, final Output output) {
        narlSerializer.write(kryo, output, novelAdjacencyAndInferredAltHaptype);
        output.writeInt(alignmentEvidence.size());
        alignmentEvidence.forEach(ca -> caSerializer.write(kryo, output, ca));
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SimpleNovelAdjacency> {
        @Override
        public void write(final Kryo kryo, final Output output, final SimpleNovelAdjacency novelAdjacencyReferenceLocations ) {
            novelAdjacencyReferenceLocations.serialize(kryo, output);
        }

        @Override
        public SimpleNovelAdjacency read(final Kryo kryo, final Input input, final Class<SimpleNovelAdjacency> klass ) {
            return new SimpleNovelAdjacency(kryo, input);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SimpleNovelAdjacency that = (SimpleNovelAdjacency) o;

        if (!novelAdjacencyAndInferredAltHaptype.equals(that.novelAdjacencyAndInferredAltHaptype)) return false;
        return alignmentEvidence.equals(that.alignmentEvidence);
    }

    @Override
    public int hashCode() {
        int result = novelAdjacencyAndInferredAltHaptype.hashCode();
        result = 31 * result + alignmentEvidence.hashCode();
        return result;
    }

    @DefaultSerializer(NovelAdjacencyAndInferredAltHaptype.Serializer.class)
    public static final class NovelAdjacencyAndInferredAltHaptype {

        private static final NovelAdjacencyReferenceLocations.Serializer localSerializer = new NovelAdjacencyReferenceLocations.Serializer();
        private final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations;
        private final byte[] altHaplotypeSequence;

        NovelAdjacencyAndInferredAltHaptype(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations,
                                            final byte[] altHaplotypeSequence) {
            this.novelAdjacencyReferenceLocations = novelAdjacencyReferenceLocations;
            this.altHaplotypeSequence = altHaplotypeSequence;
        }

        private NovelAdjacencyAndInferredAltHaptype(final Kryo kryo, final Input input) {
            novelAdjacencyReferenceLocations = localSerializer.read(kryo, input, NovelAdjacencyReferenceLocations.class);
            final boolean altSeqIsNull = input.readBoolean();
            if (altSeqIsNull) {
                altHaplotypeSequence = null;
            } else {
                final int arraySize = input.readInt();
                altHaplotypeSequence = new byte[arraySize];
                for (int i = 0 ; i < arraySize; ++i) {
                    altHaplotypeSequence[i] = input.readByte();
                }
            }
        }

        private void serialize(final Kryo kryo, final Output output) {
            localSerializer.write(kryo, output, novelAdjacencyReferenceLocations);
            if (altHaplotypeSequence==null) {
                output.writeBoolean(true);
            } else {
                output.writeBoolean(false);
                output.writeInt(altHaplotypeSequence.length);
                for (final byte b : altHaplotypeSequence) {
                    output.writeByte(b);
                }
            }
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<NovelAdjacencyAndInferredAltHaptype> {
            @Override
            public void write(final Kryo kryo, final Output output, final NovelAdjacencyAndInferredAltHaptype novelAdjacencyReferenceLocations ) {
                novelAdjacencyReferenceLocations.serialize(kryo, output);
            }

            @Override
            public NovelAdjacencyAndInferredAltHaptype read(final Kryo kryo, final Input input, final Class<NovelAdjacencyAndInferredAltHaptype> klass ) {
                return new NovelAdjacencyAndInferredAltHaptype(kryo, input);
            }
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            NovelAdjacencyAndInferredAltHaptype that = (NovelAdjacencyAndInferredAltHaptype) o;

            if (!novelAdjacencyReferenceLocations.equals(that.novelAdjacencyReferenceLocations)) return false;
            return Arrays.equals(altHaplotypeSequence, that.altHaplotypeSequence);
        }

        @Override
        public int hashCode() {
            int result = novelAdjacencyReferenceLocations.hashCode();
            result = 31 * result + Arrays.hashCode(altHaplotypeSequence);
            return result;
        }
    }
}
