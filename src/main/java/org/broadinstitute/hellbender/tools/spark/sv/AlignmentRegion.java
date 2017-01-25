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

// TODO: 11/26/16 this class is indicating there are two coordinate systems, one on the reference, the other on the assembled contig,
//       the only asymmetry should be that the sequence have exact match to the contig, but has, in theory, in-exact alignment to the reference

/**
 * A wrapper that wraps around the {@link AlnRgn} type returned by jBWA to represent one (out of the potentially multiple)
 * alignment record(s) of an locally-assembled contig.
 */
@DefaultSerializer(AlignmentRegion.Serializer.class)
class AlignmentRegion {

    static final String STRING_REP_SEPARATOR= "\t";
    static final String DUMMY_ASM_ID = "ASSEMBLY";

    final String assemblyId;
    final String contigId;

    // ref alignment info
    final SimpleInterval referenceInterval;
    final Cigar forwardStrandCigar;
    final boolean forwardStrand;
    final int mapQual;
    final int mismatches;

    // contig "exact-match" alignment info
    final int assembledContigLength;
    final int startInAssembledContig;   // 1-based, inclusive
    final int endInAssembledContig;     // 1-based, inclusive

    // TODO: the referenceInterval initialization/assignment needs +1 because of a bug in the jBWA this code depends on. fix it when Ted's binding is in.
    @VisibleForTesting
    public AlignmentRegion(final String assemblyId, final String contigId, final AlnRgn alnRgn) {
        this.assemblyId = assemblyId;
        this.contigId = contigId;

        this.forwardStrand = alnRgn.getStrand() == '+';
        final Cigar alignmentCigar = TextCigarCodec.decode(alnRgn.getCigar());
        this.forwardStrandCigar = forwardStrand ? alignmentCigar : CigarUtils.invertCigar(alignmentCigar);
        this.referenceInterval = new SimpleInterval(alnRgn.getChrom(), (int) alnRgn.getPos() + 1, (int) (alnRgn.getPos() + forwardStrandCigar.getReferenceLength()));
        this.mapQual = alnRgn.getMQual();
        this.mismatches = alnRgn.getNm();

        this.assembledContigLength = forwardStrandCigar.getReadLength() + SVVariantCallerUtils.getTotalHardClipping(forwardStrandCigar);
        this.startInAssembledContig = startOfAlignmentInContig();
        this.endInAssembledContig = endOfAlignmentInContig();
    }

    @VisibleForTesting
    public AlignmentRegion(final String assemblyId, final String contigId, final SimpleInterval referenceInterval,
                           final Cigar forwardStrandCigar, final boolean forwardStrand, final int mapQual, final int mismatches,
                           final int startInAssembledContig, final int endInAssembledContig) {
        this.assemblyId = assemblyId;
        this.contigId = contigId;
        this.referenceInterval = referenceInterval;
        this.forwardStrandCigar = forwardStrandCigar;
        this.forwardStrand = forwardStrand;
        this.mapQual = mapQual;
        this.mismatches = mismatches;
        this.assembledContigLength = forwardStrandCigar.getReadLength() + SVVariantCallerUtils.getTotalHardClipping(forwardStrandCigar);
        this.startInAssembledContig = startInAssembledContig;
        this.endInAssembledContig = endInAssembledContig;
    }

    public AlignmentRegion(final GATKRead read) {
        this.assemblyId = DUMMY_ASM_ID;
        this.contigId = read.getName();
        this.referenceInterval = new SimpleInterval(read);
        this.forwardStrand = ! read.isReverseStrand();
        this.forwardStrandCigar = forwardStrand ? read.getCigar() : CigarUtils.invertCigar(read.getCigar());
        this.mapQual = read.getMappingQuality();
        this.mismatches = read.hasAttribute("NM") ? read.getAttributeAsInteger("NM") : SVConstants.CallingStepConstants.MISSING_NM;
        this.assembledContigLength = forwardStrandCigar.getReadLength() + SVVariantCallerUtils.getTotalHardClipping(forwardStrandCigar);
        this.startInAssembledContig = startOfAlignmentInContig();
        this.endInAssembledContig = endOfAlignmentInContig();
    }

    public AlignmentRegion(final Kryo kryo, final Input input) {
        this.assemblyId = input.readString();
        this.contigId = input.readString();
        this.referenceInterval = kryo.readObject(input, SimpleInterval.class);
        this.forwardStrandCigar = TextCigarCodec.decode(input.readString());
        this.forwardStrand = input.readBoolean();
        this.mapQual = input.readInt();
        this.mismatches = input.readInt();
        this.assembledContigLength = input.readInt();
        this.startInAssembledContig = input.readInt();
        this.endInAssembledContig = input.readInt();
    }

    @VisibleForTesting
    int startOfAlignmentInContig() {
        return SVVariantCallerUtils.getNumClippedBases(true, forwardStrandCigar) + 1;
    }

    @VisibleForTesting
    int endOfAlignmentInContig() {
        return assembledContigLength - SVVariantCallerUtils.getNumClippedBases(false, forwardStrandCigar);
    }

    // TODO: 11/27/16 test
    /**
     * @return  A packed String representation of this AlignmentRegion, noticeably the fields are not separated by tab.
     *          Note that the format is NOT the same as that used in {@link #toString()}.
     */
    String toPackedString() {
        return assemblyId + "-" + contigId + ":" + startInAssembledContig + "-" + endInAssembledContig + ":" + referenceInterval.getContig() + ',' + referenceInterval.getStart() + ',' + (forwardStrand ? '+' : '-') + ',' + TextCigarCodec.encode(forwardStrandCigar) + ',' + mapQual + ',' + mismatches;
    }

    @Override
    public String toString() {
        return assemblyId +
                STRING_REP_SEPARATOR +
                contigId +
                STRING_REP_SEPARATOR +
                referenceInterval.getContig() +
                STRING_REP_SEPARATOR +
                referenceInterval.getStart() +
                STRING_REP_SEPARATOR +
                referenceInterval.getEnd() +
                STRING_REP_SEPARATOR +
                (forwardStrand ? "+" : "-") +
                STRING_REP_SEPARATOR +
                forwardStrandCigar.toString() +
                STRING_REP_SEPARATOR +
                mapQual +
                STRING_REP_SEPARATOR +
                startInAssembledContig +
                STRING_REP_SEPARATOR +
                endInAssembledContig +
                STRING_REP_SEPARATOR +
                mismatches;
    }

    /**
     * Parses fields in the same order as they were output in {@link #toString()} :
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
        final String asmId = fields[0];
        final String contigId = fields[1].replace(">", "").split(" ")[0];
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
        return new AlignmentRegion(asmId, contigId, refInterval, cigar, refStrand, mqual, mismatches, contigStart, contigEnd);
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
                Objects.equals(assemblyId, that.assemblyId) &&
                Objects.equals(contigId, that.contigId) &&
                Objects.equals(forwardStrandCigar, that.forwardStrandCigar) &&
                Objects.equals(referenceInterval, that.referenceInterval);
    }

    @Override
    public int hashCode() {
        return Objects.hash(assemblyId, contigId, forwardStrandCigar, forwardStrand, referenceInterval, mapQual, startInAssembledContig, endInAssembledContig, assembledContigLength, mismatches);
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
        output.writeString(assemblyId);
        output.writeString(contigId);
        kryo.writeObject(output, referenceInterval);
        output.writeString(TextCigarCodec.encode(forwardStrandCigar));
        output.writeBoolean(forwardStrand);
        output.writeInt(mapQual);
        output.writeInt(mismatches);
        output.writeInt(assembledContigLength);
        output.writeInt(startInAssembledContig);
        output.writeInt(endInAssembledContig);
    }
}
