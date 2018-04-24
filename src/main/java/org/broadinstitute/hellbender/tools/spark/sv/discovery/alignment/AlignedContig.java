package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import scala.Tuple2;

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

    private final String contigName;
    private final byte[] contigSequence;
    private final List<AlignmentInterval> alignmentIntervals;

    // throws if alignment interval is null
    public AlignedContig(final String contigName, final byte[] contigSequence, final List<AlignmentInterval> alignmentIntervals) {
        if (alignmentIntervals == null) {
            throw new IllegalArgumentException("AlignedContig being constructed with null alignments: " + contigName);
        }
        this.contigName = contigName;
        this.contigSequence = contigSequence;
        this.alignmentIntervals = alignmentIntervals.stream().sorted(getAlignmentIntervalComparator()).collect(Collectors.toList());
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
    }

    static Comparator<AlignmentInterval> getAlignmentIntervalComparator() {
        Comparator<AlignmentInterval> comparePos = Comparator.comparingInt(aln -> aln.startInAssembledContig);
        Comparator<AlignmentInterval> compareRefTig = Comparator.comparing(aln -> aln.referenceSpan.getContig());
        Comparator<AlignmentInterval> compareRefSpanStart = Comparator.comparingInt(aln -> aln.referenceSpan.getStart());
        return comparePos.thenComparing(compareRefTig).thenComparing(compareRefSpanStart);
    }

    boolean hasOnly2Alignments() {
        return alignmentIntervals.size() == 2;
    }

    /**
     * @return first alignment of the contig
     */
    public AlignmentInterval getHeadAlignment() {
        return alignmentIntervals.get(0);
    }

    /**
     * @return last alignment of the contig
     */
    public AlignmentInterval getTailAlignment() {
        return alignmentIntervals.get(alignmentIntervals.size() - 1);
    }

    public String getContigName() {
        return contigName;
    }

    public byte[] getContigSequence() {
        return contigSequence;
    }

    public boolean isUnmapped() {
        return alignmentIntervals.isEmpty();
    }

    public List<AlignmentInterval> getAlignments() {
        return alignmentIntervals;
    }

    @Override
    public String toString() {
        return formatContigInfo(
                new Tuple2<>(contigName, alignmentIntervals.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList())));
    }

    /**
     * Format provided {@code tigNameAndMappings} for debugging.
     */
    public static String formatContigInfo(final Tuple2<String, List<String>> tigNameAndMappings) {
        return "(" + tigNameAndMappings._1 + ", " + tigNameAndMappings._2 + ")";
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final AlignedContig that = (AlignedContig) o;

        if (!contigName.equals(that.contigName)) return false;
        if (!Arrays.equals(contigSequence, that.contigSequence)) return false;
        return alignmentIntervals.equals(that.alignmentIntervals);
    }

    @Override
    public int hashCode() {
        int result = contigName.hashCode();
        result = 31 * result + Arrays.hashCode(contigSequence);
        result = 31 * result + alignmentIntervals.hashCode();
        return result;
    }

    void serialize(final Kryo kryo, final Output output) {

        output.writeString(contigName);

        output.writeInt(contigSequence.length);
        for (final byte base : contigSequence) {
            output.writeByte(base);
        }

        output.writeInt(alignmentIntervals.size());
        alignmentIntervals.forEach(it -> it.serialize(kryo, output));

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
}
