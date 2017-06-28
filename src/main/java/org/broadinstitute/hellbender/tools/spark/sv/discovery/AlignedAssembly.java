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
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Holding necessary information about a local assembly for use in SV discovery.
 */
@DefaultSerializer(AlignedAssembly.Serializer.class)
public final class AlignedAssembly {

    public final int assemblyId;

    public final List<AlignedContig> alignedContigs;

    /**
     * Each assembled contig should have at least one such accompanying structure, or 0 when it is unmapped.
     */
    @DefaultSerializer(AlignmentInterval.Serializer.class)
    static public final class AlignmentInterval {

        public final SimpleInterval referenceInterval;
        public final int startInAssembledContig;   // 1-based, inclusive
        public final int endInAssembledContig;     // 1-based, inclusive

        public final Cigar cigarAlong5to3DirectionOfContig;

        public final boolean forwardStrand;
        public final int mapQual;
        public final int mismatches;

        @VisibleForTesting
        public AlignmentInterval(final SAMRecord samRecord) {

            final boolean isMappedReverse = samRecord.getReadNegativeStrandFlag();
            this.referenceInterval = new SimpleInterval(samRecord);
            this.startInAssembledContig = getAlignmentStartInOriginalContig(samRecord);
            this.endInAssembledContig = getAlignmentEndInOriginalContig(samRecord);

            this.cigarAlong5to3DirectionOfContig = isMappedReverse ? CigarUtils.invertCigar(samRecord.getCigar()) : samRecord.getCigar();
            this.forwardStrand = !isMappedReverse;
            this.mapQual = samRecord.getMappingQuality();
            final Integer nMismatches = samRecord.getIntegerAttribute("NM");
            this.mismatches = nMismatches==null ? SVConstants.DiscoveryStepConstants.MISSING_NM : nMismatches;
        }

        @VisibleForTesting
        public AlignmentInterval(final BwaMemAlignment alignment, final List<String> refNames, final int unclippedContigLength) {

            this.referenceInterval = new SimpleInterval(refNames.get(alignment.getRefId()), alignment.getRefStart()+1, alignment.getRefEnd()); // +1 because the BwaMemAlignment class has 0-based coordinate system
            this.forwardStrand = (alignment.getSamFlag()& SAMFlag.READ_REVERSE_STRAND.intValue())==0;
            this.cigarAlong5to3DirectionOfContig = forwardStrand ? TextCigarCodec.decode(alignment.getCigar()) : CigarUtils.invertCigar(TextCigarCodec.decode(alignment.getCigar()));
            Utils.validateArg(cigarAlong5to3DirectionOfContig.getReadLength() + SVVariantDiscoveryUtils.getTotalHardClipping(cigarAlong5to3DirectionOfContig)
                    == unclippedContigLength, "contig length provided in constructor and inferred length by computation are different: " + unclippedContigLength + "\t" + alignment.toString());

            this.mapQual = Math.max(SAMRecord.NO_MAPPING_QUALITY, alignment.getMapQual()); // BwaMemAlignment has negative mapQ for unmapped sequences, not the same as its SAMRecord conversion (see BwaMemAlignmentUtils.applyAlignment())
            this.mismatches = alignment.getNMismatches();
            if ( forwardStrand ) {
                this.startInAssembledContig = alignment.getSeqStart() + 1;
                this.endInAssembledContig = alignment.getSeqEnd();
            } else {
                this.startInAssembledContig = unclippedContigLength - alignment.getSeqEnd() + 1;
                this.endInAssembledContig = unclippedContigLength - alignment.getSeqStart();
            }
        }

        public AlignmentInterval(final SimpleInterval referenceInterval, final int startInAssembledContig, final int endInAssembledContig,
                                 final Cigar cigarAlong5to3DirectionOfContig, final boolean forwardStrand, final int mapQual, final int mismatches) {
            this.referenceInterval = referenceInterval;
            this.startInAssembledContig = startInAssembledContig;
            this.endInAssembledContig = endInAssembledContig;

            this.cigarAlong5to3DirectionOfContig = cigarAlong5to3DirectionOfContig;

            this.forwardStrand = forwardStrand;
            this.mapQual = mapQual;
            this.mismatches = mismatches;
        }

        static int getAlignmentStartInOriginalContig(final SAMRecord samRecord) {
            return SVVariantDiscoveryUtils.getNumClippedBases(!samRecord.getReadNegativeStrandFlag(), samRecord.getCigar()) + 1;
        }

        static int getAlignmentEndInOriginalContig(final SAMRecord samRecord) {
            final Cigar cigar = samRecord.getCigar();
            return cigar.getReadLength() + SVVariantDiscoveryUtils.getTotalHardClipping(cigar) - SVVariantDiscoveryUtils.getNumClippedBases(samRecord.getReadNegativeStrandFlag(), cigar);
        }

        AlignmentInterval(final Kryo kryo, final Input input) {
            final String chr = input.readString();
            final int refStart = input.readInt(),
                    refEnd = input.readInt();
            referenceInterval = new SimpleInterval(chr, refStart, refEnd);
            startInAssembledContig = input.readInt();
            endInAssembledContig = input.readInt();
            cigarAlong5to3DirectionOfContig = TextCigarCodec.decode(input.readString());
            forwardStrand = input.readBoolean();
            mapQual = input.readInt();
            mismatches = input.readInt();
        }

        void serialize(final Kryo kryo, final Output output) {
            output.writeString(referenceInterval.getContig());
            output.writeInt(referenceInterval.getStart());
            output.writeInt(referenceInterval.getEnd());
            output.writeInt(startInAssembledContig);
            output.writeInt(endInAssembledContig);
            output.writeString(TextCigarCodec.encode(cigarAlong5to3DirectionOfContig));
            output.writeBoolean(forwardStrand);
            output.writeInt(mapQual);
            output.writeInt(mismatches);
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
         * @return  A packed String representation of this alignment interval; intended for debugging or annotation usage (both requires compactified message).
         */
        String toPackedString() {
            return String.join(PACKED_STRING_REP_SEPARATOR, String.valueOf(startInAssembledContig), String.valueOf(endInAssembledContig),
                    referenceInterval.toString(), (forwardStrand ? "+" : "-"),
                    TextCigarCodec.encode(cigarAlong5to3DirectionOfContig), String.valueOf(mapQual), String.valueOf(mismatches));
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            AlignmentInterval that = (AlignmentInterval) o;

            if (startInAssembledContig != that.startInAssembledContig) return false;
            if (endInAssembledContig != that.endInAssembledContig) return false;
            if (forwardStrand != that.forwardStrand) return false;
            if (mapQual != that.mapQual) return false;
            if (mismatches != that.mismatches) return false;
            if (!referenceInterval.equals(that.referenceInterval)) return false;
            return cigarAlong5to3DirectionOfContig.equals(that.cigarAlong5to3DirectionOfContig);
        }

        @Override
        public int hashCode() {
            int result = referenceInterval.hashCode();
            result = 31 * result + startInAssembledContig;
            result = 31 * result + endInAssembledContig;
            result = 31 * result + cigarAlong5to3DirectionOfContig.hashCode();
            result = 31 * result + (forwardStrand ? 1 : 0);
            result = 31 * result + mapQual;
            result = 31 * result + mismatches;
            return result;
        }
    }

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
        for(final AlignedContig alignedContig : alignedContigs) {
            alignedContig.serialize(kryo, output);
        }
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
