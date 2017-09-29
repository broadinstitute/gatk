package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Locally assembled contig:
 *   its name
 *   its sequence as produced by the assembler (no reverse complement like in the SAM record if it maps to '-' strand), and
 *   its stripped-down alignment information.
 */
@DefaultSerializer(AlignedContig.Serializer.class)
public final class AlignedContig {

    public final String contigName;
    public final byte[] contigSequence;
    public final List<AlignmentInterval> alignmentIntervals;
    public final boolean hasEquallyGoodAlnConfigurations;

    public AlignedContig(final String contigName, final byte[] contigSequence, final List<AlignmentInterval> alignmentIntervals,
                         final boolean hasEquallyGoodAlnConfigurations) {
        this.contigName = contigName;
        this.contigSequence = contigSequence;
        this.alignmentIntervals = Utils.stream(alignmentIntervals)
                .sorted(getAlignmentIntervalComparator()).collect(Collectors.toList());
        this.hasEquallyGoodAlnConfigurations = hasEquallyGoodAlnConfigurations;
    }

    AlignedContig(final Kryo kryo, final Input input) {

        contigName = input.readString();

        final int nBases = input.readInt();
        contigSequence = new byte[nBases];
        for (int b = 0; b < nBases; ++b) {
            contigSequence[b] = input.readByte();
        }

        final int nAlignments = input.readInt();
        alignmentIntervals = new ArrayList<>(nAlignments);
        for (int i = 0; i < nAlignments; ++i) {
            alignmentIntervals.add(new AlignmentInterval(kryo, input));
        }

        hasEquallyGoodAlnConfigurations = input.readBoolean();
    }

    public static Comparator<AlignmentInterval> getAlignmentIntervalComparator() {
        Comparator<AlignmentInterval> comparePos = (AlignmentInterval a1, AlignmentInterval a2) -> Integer.compare(a1.startInAssembledContig, a2.startInAssembledContig);
        Comparator<AlignmentInterval> compareRefTig = (AlignmentInterval a1, AlignmentInterval a2) -> a1.referenceSpan.getContig().compareTo(a2.referenceSpan.getContig());
        Comparator<AlignmentInterval> compareRefSpanStart = (AlignmentInterval a1, AlignmentInterval a2) -> a1.referenceSpan.getStart() - a2.referenceSpan.getStart();
        return comparePos.thenComparing(compareRefTig).thenComparing(compareRefSpanStart);
    }
    
    void serialize(final Kryo kryo, final Output output) {

        output.writeString(contigName);

        output.writeInt(contigSequence.length);
        for (final byte base : contigSequence) {
            output.writeByte(base);
        }

        output.writeInt(alignmentIntervals.size());
        alignmentIntervals.forEach(it -> it.serialize(kryo, output));

        output.writeBoolean(hasEquallyGoodAlnConfigurations);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignedContig> {
        @Override
        public void write(final Kryo kryo, final Output output, final AlignedContig alignedContig) {
            alignedContig.serialize(kryo, output);
        }

        @Override
        public AlignedContig read(final Kryo kryo, final Input input, final Class<AlignedContig> clazz) {
            return new AlignedContig(kryo, input);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AlignedContig that = (AlignedContig) o;

        if (hasEquallyGoodAlnConfigurations != that.hasEquallyGoodAlnConfigurations) return false;
        if (!contigName.equals(that.contigName)) return false;
        if (!Arrays.equals(contigSequence, that.contigSequence)) return false;
        return alignmentIntervals.equals(that.alignmentIntervals);
    }

    @Override
    public int hashCode() {
        int result = contigName.hashCode();
        result = 31 * result + Arrays.hashCode(contigSequence);
        result = 31 * result + alignmentIntervals.hashCode();
        result = 31 * result + (hasEquallyGoodAlnConfigurations ? 1 : 0);
        return result;
    }
}
