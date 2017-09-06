package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.MISSING_NM;

/**
 * Each assembled contig should have at least one such accompanying structure, or 0 when it is unmapped.
 */
@DefaultSerializer(AlignmentInterval.Serializer.class)
public final class AlignmentInterval {

    public final SimpleInterval referenceSpan;
    public final int startInAssembledContig;   // 1-based, inclusive
    public final int endInAssembledContig;     // 1-based, inclusive

    public final Cigar cigarAlong5to3DirectionOfContig;

    public final boolean forwardStrand;
    public final int mapQual;
    public final int mismatches;
    public final int alnScore;

    // if any of the following boolean fields are true, fields "mapQual", "mismatches", "alnScore" should be viewed
    // with care as they were simply copied from the original alignment (not for "mismatches"),
    // which after the split are wrong (we didn't recompute them because that would require expensive SW re-alignment)
    public final boolean isFromSplitGapAlignment;
    public final boolean hasUndergoneOverlapRemoval;

    @VisibleForTesting
    public AlignmentInterval(final SAMRecord samRecord) {

        final boolean isMappedReverse = samRecord.getReadNegativeStrandFlag();
        this.referenceSpan = new SimpleInterval(samRecord);
        this.startInAssembledContig = getAlignmentStartInOriginalContig(samRecord);
        this.endInAssembledContig = getAlignmentEndInOriginalContig(samRecord);

        this.cigarAlong5to3DirectionOfContig = isMappedReverse ? CigarUtils.invertCigar(samRecord.getCigar())
                                                               : samRecord.getCigar();
        this.forwardStrand = !isMappedReverse;
        this.mapQual = samRecord.getMappingQuality();
        final Integer nMismatches = samRecord.getIntegerAttribute("NM");
        this.mismatches = nMismatches==null ? MISSING_NM : nMismatches;

        int alignerScore;
        try {
            alignerScore = samRecord.getIntegerAttribute("AS");
        } catch (final RuntimeException ex) {
            alignerScore = -1;
        }
        this.alnScore = alignerScore;
        this.isFromSplitGapAlignment = false;
        this.hasUndergoneOverlapRemoval = false;
    }

    @VisibleForTesting
    public AlignmentInterval(final BwaMemAlignment alignment, final List<String> refNames, final int unclippedContigLength) {

        // +1 because the BwaMemAlignment class has 0-based coordinate system
        this.referenceSpan = new SimpleInterval(refNames.get(alignment.getRefId()),
                                            alignment.getRefStart() + 1, alignment.getRefEnd());
        this.forwardStrand = 0==(alignment.getSamFlag() & SAMFlag.READ_REVERSE_STRAND.intValue());
        this.cigarAlong5to3DirectionOfContig = forwardStrand ? TextCigarCodec.decode(alignment.getCigar())
                                                               : CigarUtils.invertCigar(TextCigarCodec.decode(alignment.getCigar()));
        Utils.validateArg(
                cigarAlong5to3DirectionOfContig.getReadLength() + SvCigarUtils.getTotalHardClipping(cigarAlong5to3DirectionOfContig)
                        == unclippedContigLength,
                "contig length provided in constructor and inferred length by computation are different: " +
                        unclippedContigLength + "\t" + alignment.toString());

        // BwaMemAlignment has negative mapQ for unmapped sequences, not the same as its SAMRecord conversion
        // (see BwaMemAlignmentUtils.applyAlignment())
        this.mapQual = Math.max(SAMRecord.NO_MAPPING_QUALITY, alignment.getMapQual());
        this.mismatches = alignment.getNMismatches();
        if (forwardStrand) {
            this.startInAssembledContig = alignment.getSeqStart() + 1;
            this.endInAssembledContig = alignment.getSeqEnd();
        } else {
            this.startInAssembledContig = unclippedContigLength - alignment.getSeqEnd() + 1;
            this.endInAssembledContig = unclippedContigLength - alignment.getSeqStart();
        }

        this.alnScore = alignment.getAlignerScore();
        this.isFromSplitGapAlignment = false;
        this.hasUndergoneOverlapRemoval = false;
    }

    public AlignmentInterval(final SimpleInterval referenceSpan, final int startInAssembledContig, final int endInAssembledContig,
                             final Cigar cigarAlong5to3DirectionOfContig, final boolean forwardStrand,
                             final int mapQual, final int mismatches, final int alignerScore,
                             final boolean isFromSplitGapAlignment, final boolean hasUndergoneOverlapRemoval) {
        this.referenceSpan = referenceSpan;
        this.startInAssembledContig = startInAssembledContig;
        this.endInAssembledContig = endInAssembledContig;

        this.cigarAlong5to3DirectionOfContig = cigarAlong5to3DirectionOfContig;

        this.forwardStrand = forwardStrand;
        this.mapQual = mapQual;
        this.mismatches = mismatches;
        this.alnScore = alignerScore;
        this.isFromSplitGapAlignment = isFromSplitGapAlignment;
        this.hasUndergoneOverlapRemoval = hasUndergoneOverlapRemoval;
    }

    /**
     * @return the number of bases of overlap between two alignment regions overlap on the locally-assembled contig they originate from.
     *          Mostly useful for computing micro-homologyForwardStrandRep.
     */
    @VisibleForTesting
    public static int overlapOnContig(final AlignmentInterval one, final AlignmentInterval two) {
        return Math.max(0,
                Math.min(one.endInAssembledContig + 1, two.endInAssembledContig + 1)
                        - Math.max(one.startInAssembledContig, two.startInAssembledContig)
        );
    }

    /**
     * Computes overlap between reference span of the two input alignment intervals.
     */
    static int overlapOnRefSpan(final AlignmentInterval one, final AlignmentInterval two) {

        if ( !one.referenceSpan.getContig().equals(two.referenceSpan.getContig()) ) return  0;

        // dummy number for chr to be used in constructing SVInterval, since 2 input AI's both map to the same chr by this point
        final int dummyChr = -1;
        final SVInterval intOne = new SVInterval(dummyChr, one.referenceSpan.getStart(), one.referenceSpan.getEnd() + 1),
                         intTwo = new SVInterval(dummyChr, two.referenceSpan.getStart(), two.referenceSpan.getEnd() + 1);

        return intOne.overlapLen(intTwo);
    }

    static int getAlignmentStartInOriginalContig(final SAMRecord samRecord) {
        return SvCigarUtils.getNumClippedBases(!samRecord.getReadNegativeStrandFlag(), samRecord.getCigar()) + 1;
    }

    static int getAlignmentEndInOriginalContig(final SAMRecord samRecord) {
        final Cigar cigar = samRecord.getCigar();
        return cigar.getReadLength() + SvCigarUtils.getTotalHardClipping(cigar) -
                SvCigarUtils.getNumClippedBases(samRecord.getReadNegativeStrandFlag(), cigar);
    }

    AlignmentInterval(final Kryo kryo, final Input input) {
        final String chr   = input.readString();
        final int refStart = input.readInt(),
                  refEnd   = input.readInt();
        referenceSpan = new SimpleInterval(chr, refStart, refEnd);
        startInAssembledContig = input.readInt();
        endInAssembledContig = input.readInt();
        cigarAlong5to3DirectionOfContig = TextCigarCodec.decode(input.readString());
        forwardStrand = input.readBoolean();
        mapQual = input.readInt();
        mismatches = input.readInt();
        alnScore = input.readInt();
        isFromSplitGapAlignment = input.readBoolean();
        hasUndergoneOverlapRemoval = input.readBoolean();
    }

    void serialize(final Kryo kryo, final Output output) {
        output.writeString(referenceSpan.getContig());
        output.writeInt(referenceSpan.getStart());
        output.writeInt(referenceSpan.getEnd());
        output.writeInt(startInAssembledContig);
        output.writeInt(endInAssembledContig);
        output.writeString(TextCigarCodec.encode(cigarAlong5to3DirectionOfContig));
        output.writeBoolean(forwardStrand);
        output.writeInt(mapQual);
        output.writeInt(mismatches);
        output.writeInt(alnScore);
        output.writeBoolean(isFromSplitGapAlignment);
        output.writeBoolean(hasUndergoneOverlapRemoval);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignmentInterval> {
        @Override
        public void write( final Kryo kryo, final Output output, final AlignmentInterval alignmentInterval){
            alignmentInterval.serialize(kryo, output);
        }

        @Override
        public AlignmentInterval read(final Kryo kryo, final Input input, final Class<AlignmentInterval> clazz ) {
            return new AlignmentInterval(kryo, input);
        }
    }

    static final String PACKED_STRING_REP_SEPARATOR = "_";
    /**
     * @return  A packed String representation of this alignment interval; intended for debugging or annotation usage
     *          (both requires compactified message).
     */
    public String toPackedString() {
        return String.join(PACKED_STRING_REP_SEPARATOR, String.valueOf(startInAssembledContig),
                String.valueOf(endInAssembledContig), referenceSpan.toString(), (forwardStrand ? "+" : "-"),
                TextCigarCodec.encode(cigarAlong5to3DirectionOfContig),
                String.valueOf(mapQual), String.valueOf(mismatches), String.valueOf(alnScore),
                (isFromSplitGapAlignment ? "s" : "o"), (hasUndergoneOverlapRemoval ? "h" : "nh"));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AlignmentInterval interval = (AlignmentInterval) o;

        if (startInAssembledContig != interval.startInAssembledContig) return false;
        if (endInAssembledContig != interval.endInAssembledContig) return false;
        if (forwardStrand != interval.forwardStrand) return false;
        if (mapQual != interval.mapQual) return false;
        if (mismatches != interval.mismatches) return false;
        if (alnScore != interval.alnScore) return false;
        if (isFromSplitGapAlignment != interval.isFromSplitGapAlignment) return false;
        if (hasUndergoneOverlapRemoval != interval.hasUndergoneOverlapRemoval) return false;
        if (!referenceSpan.equals(interval.referenceSpan)) return false;
        return cigarAlong5to3DirectionOfContig.equals(interval.cigarAlong5to3DirectionOfContig);
    }

    @Override
    public int hashCode() {
        int result = referenceSpan.hashCode();
        result = 31 * result + startInAssembledContig;
        result = 31 * result + endInAssembledContig;
        result = 31 * result + cigarAlong5to3DirectionOfContig.hashCode();
        result = 31 * result + (forwardStrand ? 1 : 0);
        result = 31 * result + mapQual;
        result = 31 * result + mismatches;
        result = 31 * result + alnScore;
        result = 31 * result + (isFromSplitGapAlignment ? 1 : 0);
        result = 31 * result + (hasUndergoneOverlapRemoval ? 1 : 0);
        return result;
    }
}
