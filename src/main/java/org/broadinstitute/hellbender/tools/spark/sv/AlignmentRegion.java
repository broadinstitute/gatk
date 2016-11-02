package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.github.lindenb.jbwa.jni.AlnRgn;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Objects;

/**
 * A wrapper that wraps around the {@link AlnRgn} type returned by jBWA to represent one (out of the potentially multiple)
 * alignment record(s) of an locally-assembled contig.
 */
@DefaultSerializer(AlignmentRegion.Serializer.class)
class AlignmentRegion {

    final String contigId;
    final String assemblyId;
    final Cigar forwardStrandCigar;
    final boolean forwardStrand;
    final SimpleInterval referenceInterval;
    final int mapQual;
    final int startInAssembledContig;   // 1-based, inclusive
    final int endInAssembledContig;     // 1-based, inclusive
    final int assembledContigLength;
    final int mismatches;

    // TODO: the referenceInterval initialization/assignment needs +1 because of a bug in the jBWA this code depends on. fix it when Ted's binding is in.
    public AlignmentRegion(final String assemblyId, final String contigId, final AlnRgn alnRgn) {
        this.contigId = contigId;
        this.assemblyId = assemblyId;
        this.forwardStrand = alnRgn.getStrand() == '+';
        final Cigar alignmentCigar = TextCigarCodec.decode(alnRgn.getCigar());
        this.forwardStrandCigar = forwardStrand ? alignmentCigar : CigarUtils.invertCigar(alignmentCigar);
        this.referenceInterval = new SimpleInterval(alnRgn.getChrom(), (int) alnRgn.getPos() + 1, (int) (alnRgn.getPos() + forwardStrandCigar.getReferenceLength()));
        this.mapQual = alnRgn.getMQual();
        this.assembledContigLength = forwardStrandCigar.getReadLength() + SVVariantCallerUtils.getTotalHardClipping(forwardStrandCigar);
        this.startInAssembledContig = startOfAlignmentInContig();
        this.endInAssembledContig = endOfAlignmentInContig();
        this.mismatches = alnRgn.getNm();
    }

    public AlignmentRegion(final String assemblyId, final String contigId, final Cigar forwardStrandCigar, final boolean forwardStrand, final SimpleInterval referenceInterval, final int mapQual, final int startInAssembledContig, final int endInAssembledContig, final int mismatches) {
        this.contigId = contigId;
        this.assemblyId = assemblyId;
        this.forwardStrandCigar = forwardStrandCigar;
        this.forwardStrand = forwardStrand;
        this.referenceInterval = referenceInterval;
        this.mapQual = mapQual;
        this.startInAssembledContig = startInAssembledContig;
        this.endInAssembledContig = endInAssembledContig;
        this.assembledContigLength = forwardStrandCigar.getReadLength() + SVVariantCallerUtils.getTotalHardClipping(forwardStrandCigar);
        this.mismatches = mismatches;
    }

    public AlignmentRegion(final GATKRead read) {
        this.assemblyId = null;
        this.contigId = read.getName();
        this.forwardStrand = ! read.isReverseStrand();
        this.forwardStrandCigar = forwardStrand ? read.getCigar() : CigarUtils.invertCigar(read.getCigar());
        this.referenceInterval = new SimpleInterval(read);
        this.assembledContigLength = forwardStrandCigar.getReadLength() + SVVariantCallerUtils.getTotalHardClipping(forwardStrandCigar);
        this.startInAssembledContig = startOfAlignmentInContig();
        this.endInAssembledContig = endOfAlignmentInContig();
        this.mapQual = read.getMappingQuality();
        this.mismatches = read.hasAttribute("NM") ? read.getAttributeAsInteger("NM") : 0;
    }

    public AlignmentRegion(final Kryo kryo, final Input input) {
        this.contigId = input.readString();
        this.assemblyId = input.readString();
        this.forwardStrandCigar = TextCigarCodec.decode(input.readString());
        this.forwardStrand = input.readBoolean();
        this.referenceInterval = kryo.readObject(input, SimpleInterval.class);
        this.mapQual = input.readInt();
        this.startInAssembledContig = input.readInt();
        this.endInAssembledContig = input.readInt();
        this.assembledContigLength = input.readInt();
        this.mismatches = input.readInt();
    }

    @VisibleForTesting
    int startOfAlignmentInContig() {
        return SVVariantCallerUtils.getNumClippedBases(true, forwardStrandCigar) + 1;
    }

    @VisibleForTesting
    int endOfAlignmentInContig() {
        return assembledContigLength - SVVariantCallerUtils.getNumClippedBases(false, forwardStrandCigar);
    }

    @Override
    public String toString() {
        return assemblyId +
                "\t" +
                contigId +
                "\t" +
                referenceInterval.getContig() +
                "\t" +
                referenceInterval.getStart() +
                "\t" +
                referenceInterval.getEnd() +
                "\t" +
                (forwardStrand ? "+" : "-") +
                "\t" +
                forwardStrandCigar.toString() +
                "\t" +
                mapQual +
                "\t" +
                startInAssembledContig +
                "\t" +
                endInAssembledContig +
                "\t" +
                mismatches;
    }

    /**
     * Parses fields in the same order as they were output in {@link ContigsCollection#toString()} :
     *
     * assemblyId
     * contigId
     * refContig
     * referenceInterval.getContig
     * referenceInterval.getStart
     * referenceInterval.getEnd
     * forwardStrand
     * forwardStrandCigar
     * startInAssembledContig
     * endInAssembledContig
     * mismatches
     */
    @VisibleForTesting
    static AlignmentRegion fromString(final String[] fields) {
        final String breakpointId = fields[0];
        final String contigId = fields[1];
        final String refContig = fields[2];
        final Integer refStart = Integer.valueOf(fields[3]);
        final Integer refEnd = Integer.valueOf(fields[4]);
        final SimpleInterval refInterval = new SimpleInterval(refContig, refStart, refEnd);
        final boolean refStrand = ("+".equals(fields[5]));
        final Cigar cigar = TextCigarCodec.decode(fields[6]);
        final int mqual = Integer.valueOf(fields[7]);
        final int contigStart = Integer.valueOf(fields[8]);
        final int contigEnd = Integer.valueOf(fields[9]);
        final int mismatches = Integer.valueOf(fields[10]);
        return new AlignmentRegion(breakpointId, contigId, cigar, refStrand, refInterval, mqual, contigStart, contigEnd, mismatches);
    }

    /**
     * input format is the text representation of an alignment region
     * @param alignedAssembledContigLine An input line with the tab-separated fields of an alignment region
     * @return A tuple with the breakpoint ID and string representation of an ChimericAlignment, or an empty iterator if the line did not have two comma-separated values
     */
    @VisibleForTesting
    static AlignmentRegion parseAlignedAssembledContigLine(final String alignedAssembledContigLine) {
        final String[] split = alignedAssembledContigLine.split("\t", -1);
        return AlignmentRegion.fromString(split);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        final AlignmentRegion that = (AlignmentRegion) o;
        return forwardStrand == that.forwardStrand &&
                mapQual == that.mapQual &&
                startInAssembledContig == that.startInAssembledContig &&
                endInAssembledContig == that.endInAssembledContig &&
                assembledContigLength == that.assembledContigLength &&
                mismatches == that.mismatches &&
                Objects.equals(contigId, that.contigId) &&
                Objects.equals(assemblyId, that.assemblyId) &&
                Objects.equals(forwardStrandCigar, that.forwardStrandCigar) &&
                Objects.equals(referenceInterval, that.referenceInterval);
    }

    @Override
    public int hashCode() {
        return Objects.hash(contigId, assemblyId, forwardStrandCigar, forwardStrand, referenceInterval, mapQual, startInAssembledContig, endInAssembledContig, assembledContigLength, mismatches);
    }

    String toPackedString() {
        return assemblyId + "-" + contigId + ":" + startInAssembledContig + "-" + endInAssembledContig + ":" + referenceInterval.getContig() + ',' + referenceInterval.getStart() + ',' + (forwardStrand ? '+' : '-') + ',' + TextCigarCodec.encode(forwardStrandCigar) + ',' + mapQual + ',' + mismatches;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignmentRegion> {
        @Override
        public void write(final Kryo kryo, final Output output, final AlignmentRegion alignmentRegion ) {
            alignmentRegion.serialize(kryo, output);
        }

        @Override
        public AlignmentRegion read(final Kryo kryo, final Input input, final Class<AlignmentRegion> klass ) {
            return new AlignmentRegion(kryo, input);
        }
    }

    private void serialize(final Kryo kryo, final Output output) {
        output.writeString(contigId);
        output.writeString(assemblyId);
        output.writeString(TextCigarCodec.encode(forwardStrandCigar));
        output.writeBoolean(forwardStrand);
        kryo.writeObject(output, referenceInterval);
        output.writeInt(mapQual);
        output.writeInt(startInAssembledContig);
        output.writeInt(endInAssembledContig);
        output.writeInt(assembledContigLength);
        output.writeInt(mismatches);
    }
}
