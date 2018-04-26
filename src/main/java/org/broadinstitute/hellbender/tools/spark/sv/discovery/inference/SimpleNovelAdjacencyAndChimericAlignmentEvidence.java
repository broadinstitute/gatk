package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;

@DefaultSerializer(SimpleNovelAdjacencyAndChimericAlignmentEvidence.Serializer.class)
public final class SimpleNovelAdjacencyAndChimericAlignmentEvidence {
    private static final NovelAdjacencyAndAltHaplotype.Serializer narlSerializer = new NovelAdjacencyAndAltHaplotype.Serializer();
    private static final SimpleChimeraAndNCAMstring.Serializer alignmentEvidenceSerializer = new SimpleChimeraAndNCAMstring.Serializer();

    private final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype;
    private final List<SimpleChimeraAndNCAMstring> alignmentEvidence;

    public NovelAdjacencyAndAltHaplotype getNovelAdjacencyReferenceLocations() {
        return novelAdjacencyAndAltHaplotype;
    }
    public byte[] getAltHaplotypeSequence() {
        return novelAdjacencyAndAltHaplotype.getAltHaplotypeSequence();
    }
    public List<SimpleChimeraAndNCAMstring> getAlignmentEvidence() {
        return alignmentEvidence;
    }

    SimpleNovelAdjacencyAndChimericAlignmentEvidence(final NovelAdjacencyAndAltHaplotype novelAdjacencyReferenceLocations,
                                                     final Iterable<SimpleChimeraAndNCAMstring> alignmentEvidence) {
        this.novelAdjacencyAndAltHaplotype = Utils.nonNull( novelAdjacencyReferenceLocations );
        this.alignmentEvidence = Lists.newArrayList( Utils.nonNull(alignmentEvidence) );
    }

    private SimpleNovelAdjacencyAndChimericAlignmentEvidence(final Kryo kryo, final Input input) {
        novelAdjacencyAndAltHaplotype = narlSerializer.read(kryo, input, NovelAdjacencyAndAltHaplotype.class);
        final int evidenceCount = input.readInt();
        alignmentEvidence = new ArrayList<>(evidenceCount);
        for (int i = 0; i < evidenceCount; ++i) {
            alignmentEvidence.add( alignmentEvidenceSerializer.read(kryo, input, SimpleChimeraAndNCAMstring.class) );
        }
    }

    private void serialize(final Kryo kryo, final Output output) {
        narlSerializer.write(kryo, output, novelAdjacencyAndAltHaplotype);
        output.writeInt(alignmentEvidence.size());
        alignmentEvidence.forEach(ev -> alignmentEvidenceSerializer.write(kryo, output, ev));
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SimpleNovelAdjacencyAndChimericAlignmentEvidence> {
        @Override
        public void write(final Kryo kryo, final Output output, final SimpleNovelAdjacencyAndChimericAlignmentEvidence novelAdjacencyReferenceLocations ) {
            novelAdjacencyReferenceLocations.serialize(kryo, output);
        }

        @Override
        public SimpleNovelAdjacencyAndChimericAlignmentEvidence read(final Kryo kryo, final Input input, final Class<SimpleNovelAdjacencyAndChimericAlignmentEvidence> klass ) {
            return new SimpleNovelAdjacencyAndChimericAlignmentEvidence(kryo, input);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SimpleNovelAdjacencyAndChimericAlignmentEvidence that = (SimpleNovelAdjacencyAndChimericAlignmentEvidence) o;

        if (!novelAdjacencyAndAltHaplotype.equals(that.novelAdjacencyAndAltHaplotype)) return false;
        return alignmentEvidence.equals(that.alignmentEvidence);
    }

    @Override
    public int hashCode() {
        int result = novelAdjacencyAndAltHaplotype.hashCode();
        result = 31 * result + alignmentEvidence.hashCode();
        return result;
    }

    /**
     * A simple struct holding information of simple chimera of a particular contig and
     * its good mapping to a non-canonical chromosome, if it has one.
     */
    @DefaultSerializer(SimpleChimeraAndNCAMstring.Serializer.class)
    public static final class SimpleChimeraAndNCAMstring {
        private static final SimpleChimera.Serializer caSerializer = new SimpleChimera.Serializer();

        public final SimpleChimera simpleChimera;
        public final String goodNonCanonicalMappingSATag;

        public SimpleChimeraAndNCAMstring(final SimpleChimera simpleChimera, final String goodNonCanonicalMappingSATag) {
            this.simpleChimera = simpleChimera;
            this.goodNonCanonicalMappingSATag = goodNonCanonicalMappingSATag;
        }

        SimpleChimeraAndNCAMstring(final Kryo kryo, final Input input) {
            simpleChimera = caSerializer.read(kryo, input, SimpleChimera.class);
            goodNonCanonicalMappingSATag = input.readString();
        }

        private void serialize(final Kryo kryo, final Output output) {
            caSerializer.write(kryo, output, simpleChimera);
            output.writeString(goodNonCanonicalMappingSATag);
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final SimpleChimeraAndNCAMstring that = (SimpleChimeraAndNCAMstring) o;

            if (!simpleChimera.equals(that.simpleChimera)) return false;
            return goodNonCanonicalMappingSATag.equals(that.goodNonCanonicalMappingSATag);
        }

        @Override
        public int hashCode() {
            int result = simpleChimera.hashCode();
            result = 31 * result + goodNonCanonicalMappingSATag.hashCode();
            return result;
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SimpleChimeraAndNCAMstring> {
            @Override
            public void write(final Kryo kryo, final Output output, final SimpleChimeraAndNCAMstring alignmentEvidence) {
                alignmentEvidence.serialize(kryo, output);
            }

            @Override
            public SimpleChimeraAndNCAMstring read(final Kryo kryo, final Input input, final Class<SimpleChimeraAndNCAMstring> klass) {
                return new SimpleChimeraAndNCAMstring(kryo, input);
            }
        }
    }
}
