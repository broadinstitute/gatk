package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Each assembled contig should have at least one such accompanying structure, or 0 when it is unmapped.
 */
@DefaultSerializer(AlignmentInterval.Serializer.class)
public final class AlignmentInterval {

    /**
     * Indicates that no NM (number of mismatches) annotation is present.
     */
    public static final int NO_NM = -1;

    /**
     * Indicates that no AS (alignment score) annotation is present.
     */
    public static final int NO_AS = -1;

    public final SimpleInterval referenceSpan;
    public final int startInAssembledContig;   // 1-based, inclusive
    public final int endInAssembledContig;     // 1-based, inclusive

    public final Cigar cigarAlong5to3DirectionOfContig;

    public final boolean forwardStrand;

    /**
     * Mapping quality, {@link SAMRecord#NO_MAPPING_QUALITY} if unspecified.
     */
    public final int mapQual;

    /**
     * Number of mismatches, {@link #NO_NM} if unknown or unspecified.
     */
    public final int mismatches;

    /**
     * Alignment score, {@link #NO_AS} if unknown or unspecified.
     */
    public final int alnScore;

    // if any of the following boolean fields are true, fields "mapQual", "mismatches", "alnScore" should be viewed
    // with care as they were simply copied from the original alignment (not for "mismatches"),
    // which after the split are wrong (we didn't recompute them because that would require expensive SW re-alignment)

    public final boolean isFromSplitGapAlignment;
    public final boolean hasUndergoneOverlapRemoval;

    /**
     * Compose an alignment interval instance from a SAM supplementary alignment formatted string.
     * <p>
     *     The input string format is:
     *     <pre>
     *         contig-name,start,strand,cigar,mq,nm,as
     *     </pre>
     *     where nm (number of mismatches) and as (alignment score) might be absent. The strand
     *     is symbolized as '+' for the forward-strand and '-' for the reverse-strand.
     *     <p>An example with different degrees of completeness:</p>
     *     <pre>
     *         chr10,1241241,+,10S1313M45I14M100H,30,5,66
     *         chr10,1241241,+,10S1313M45I14M100H,30,5
     *         chr10,1241241,+,10S1313M45I14M100H,30
     *     </pre>
     *
     * </p>
     * @param samSAtagString the input string.
     * @throws IllegalArgumentException if {@code str} is {@code null} or it does not look like a valid
     *                                  SA string.
     */
    public AlignmentInterval(final String samSAtagString) {
        Utils.nonNull(samSAtagString, "input str cannot be null");
        final String[] parts = samSAtagString.replaceAll(";$", "").split(",");
        if (parts.length < 5) {
            throw new IllegalArgumentException("the input SA string at least must contain 5 parts: " + samSAtagString);
        }

        int nextPartsIndex = 0; // holds the next element index in parts to be parsed.

        final String referenceContig = parts[nextPartsIndex++];
        final int start = Integer.parseInt(parts[nextPartsIndex++]);
        final Strand strand = parseStrand(parts[nextPartsIndex++]);
        final boolean forwardStrand = strand == Strand.POSITIVE;
        final Cigar originalCigar = TextCigarCodec.decode(parts[nextPartsIndex++]);
        final Cigar cigar = forwardStrand ? originalCigar : CigarUtils.invertCigar(originalCigar);
        final int mappingQuality = ParamUtils.inRange(parseInt(parts[nextPartsIndex++]), 0, 255, "the mapping quality must be in the range [0, 255]");
        final int mismatches = parts.length > nextPartsIndex ? ParamUtils.isPositiveOrZero(parseInt(parts[nextPartsIndex++]), "the number of mismatches cannot be negative") : NO_NM;
        final int alignmentScore = parts.length > nextPartsIndex ? ParamUtils.isPositiveOrZero(parseInt(parts[nextPartsIndex]), "the alignment score cannot be negative") : NO_AS;

        this.referenceSpan = new SimpleInterval(referenceContig, start,
                Math.max(start, cigar.getReferenceLength() + start - 1));
        this.startInAssembledContig = 1 + CigarUtils.countLeftClippedBases(cigar);
        this.endInAssembledContig = CigarUtils.countUnclippedReadBases(cigar) - CigarUtils.countRightClippedBases(cigar);
        this.mapQual = mappingQuality;
        this.mismatches = mismatches;
        this.alnScore = alignmentScore;
        this.forwardStrand = forwardStrand;
        this.cigarAlong5to3DirectionOfContig = cigar;
        this.isFromSplitGapAlignment = false;
        this.hasUndergoneOverlapRemoval = false;
    }

    private static int parseInt(final String str) {
        try {
            return Integer.parseInt(str);
        } catch (NumberFormatException ex) {
            throw new IllegalArgumentException("not a valid integer in: " + str, ex);
        }
    }

    /**
     * Parses a strand text into a {@link Strand} object.
     * <p>
     * Makes sure that the result isn't the valid strand {@link Strand#NONE}.
     * </p>
     *
     * @param strandText the input strand text.
     * @return never {@code null}, either {@link Strand#POSITIVE} or {@link Strand#NEGATIVE}.
     * @throws IllegalArgumentException if {@code strandText} is not a valid strand representation
     *                                  or makes reference to {@link Strand#NONE}.
     */
    private static Strand parseStrand(final String strandText) {
        final Strand strand = Strand.toStrand(strandText);
        if (strand == Strand.NONE) {
            throw new IllegalArgumentException("the input strand cannot be: " + strandText);
        }
        return strand;
    }

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
        this.mismatches = ReadUtils.getOptionalIntAttribute(samRecord, SAMTag.NM.name()).orElse(NO_NM);
        this.alnScore = ReadUtils.getOptionalIntAttribute(samRecord, SAMTag.AS.name()).orElse(NO_AS);
        this.isFromSplitGapAlignment = false;
        this.hasUndergoneOverlapRemoval = false;
    }

    /**
     * Constructs an alignment interval that reflects on the mapping properties of a {@link GATKRead} instance.
     *
     * @param read the target read.
     * @throws IllegalArgumentException if {@code read} is {@code null}.
     */
    public AlignmentInterval(final GATKRead read) {
        Utils.nonNull(read, "the input read cannot be null");
        final boolean isMappedReverse = read.isReverseStrand();
        this.referenceSpan = new SimpleInterval(read);
        this.startInAssembledContig = ReadUtils.getFirstAlignedReadPosition(read);
        this.endInAssembledContig = ReadUtils.getLastAlignedReadPosition(read);
        this.cigarAlong5to3DirectionOfContig = isMappedReverse ? CigarUtils.invertCigar(read.getCigar()) : read.getCigar();
        this.forwardStrand = !isMappedReverse;
        this.mapQual = read.getMappingQuality();
        this.mismatches = ReadUtils.getOptionalIntAttribute(read, SAMTag.NM.name()).orElse(NO_NM);
        this.alnScore = ReadUtils.getOptionalIntAttribute(read, SAMTag.AS.name()).orElse(NO_AS);
        this.isFromSplitGapAlignment = false;
        this.hasUndergoneOverlapRemoval = false;
    }

    @VisibleForTesting
    public AlignmentInterval(final BwaMemAlignment alignment, final List<String> refNames, final int unclippedContigLength) {

        // +1 because the BwaMemAlignment class has 0-based coordinate system
        this.referenceSpan = new SimpleInterval(refNames.get(alignment.getRefId()),
                alignment.getRefStart() + 1, alignment.getRefEnd());
        this.forwardStrand = 0 == (alignment.getSamFlag() & SAMFlag.READ_REVERSE_STRAND.intValue());
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
     * Mostly useful for computing micro-homologyForwardStrandRep.
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

        if (!one.referenceSpan.getContig().equals(two.referenceSpan.getContig())) return 0;

        // dummy number for chr to be used in constructing SVInterval, since 2 input AI's both map to the same chr by this point
        final int dummyChr = -1;
        final SVInterval intOne = new SVInterval(dummyChr, one.referenceSpan.getStart(), one.referenceSpan.getEnd() + 1),
                intTwo = new SVInterval(dummyChr, two.referenceSpan.getStart(), two.referenceSpan.getEnd() + 1);

        return intOne.overlapLen(intTwo);
    }

    private static int getAlignmentStartInOriginalContig(final SAMRecord samRecord) {
        return SvCigarUtils.getNumClippedBases(!samRecord.getReadNegativeStrandFlag(), samRecord.getCigar()) + 1;
    }

    private static int getAlignmentEndInOriginalContig(final SAMRecord samRecord) {
        final Cigar cigar = samRecord.getCigar();
        return cigar.getReadLength() + SvCigarUtils.getTotalHardClipping(cigar) -
                SvCigarUtils.getNumClippedBases(samRecord.getReadNegativeStrandFlag(), cigar);
    }

    AlignmentInterval(final Kryo kryo, final Input input) {
        final String chr = input.readString();
        final int refStart = input.readInt(),
                refEnd = input.readInt();
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
        public void write(final Kryo kryo, final Output output, final AlignmentInterval alignmentInterval) {
            alignmentInterval.serialize(kryo, output);
        }

        @Override
        public AlignmentInterval read(final Kryo kryo, final Input input, final Class<AlignmentInterval> clazz) {
            return new AlignmentInterval(kryo, input);
        }
    }

    @VisibleForTesting
    static final String PACKED_STRING_REP_SEPARATOR = "_";

    /**
     * @return A packed String representation of this alignment interval; intended for debugging or annotation usage
     * (both requires compactified message).
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

    /**
     * Returns a {@link SAMRecord} instance that reflects this alignment interval given the
     * output {@link SAMFileHeader} and the enclosing sequence bases.
     *
     * @param header          the returned record header.
     * @param name            name to given to the resulting record.
     * @param unclippedBases  the enclosing contig/read bases.
     * @param hardClip        whether clippings must be hard ones ({@code true}) or soft-ones ({@code false}).
     * @param otherFlags      values for flags whose information in not contained in a {@link AlignmentInterval}.
     *                        for flags that are somehow covered by {@link AlignmentInterval} we ignored this argument values.
     * @param otherAttributes values for attributes other than the ones this instance has values for; attributes
     *                        for which this instance has value for are ignored.
     * @return never {@code null}.
     * @throws IllegalArgumentException if either {@code header} or {@code contig} is {@code null}.
     */
    public SAMRecord toSAMRecord(final SAMFileHeader header, final String name,
                                 final byte[] unclippedBases, final boolean hardClip, final int otherFlags,
                                 final Collection<? extends SAMRecord.SAMTagAndValue> otherAttributes) {
        Utils.nonNull(header, "the input header cannot be null");
        final SAMRecord result = new SAMRecord(header);
        if (hardClip && SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(otherFlags) && SAMFlag.SUPPLEMENTARY_ALIGNMENT.isUnset(otherFlags)) {
            throw new IllegalArgumentException("you cannot request hard-clipping on a primary non-supplementary alignment record");
        }

        result.setReadName(name);
        result.setFlags(otherFlags);
        result.setReadNegativeStrandFlag(!forwardStrand);

        // taking care of the bases;
        final byte[] bases = unclippedBases == null
                ? null
                : hardClip ? Arrays.copyOfRange(unclippedBases, startInAssembledContig - 1, endInAssembledContig)
                           : unclippedBases.clone();
        if (!forwardStrand && bases != null) {
            SequenceUtil.reverseComplement(bases);
        }
        result.setReadBases(bases);

        // taking care of the cigar.
        final Cigar cigar = forwardStrand ? this.cigarAlong5to3DirectionOfContig :
                CigarUtils.invertCigar(this.cigarAlong5to3DirectionOfContig);

        result.setCigar(softOrHardReclip(cigar, hardClip ? CigarOperator.H : CigarOperator.S));

        result.setReferenceName(referenceSpan.getContig());
        result.setAlignmentStart(referenceSpan.getStart());
        result.setMappingQuality(this.mapQual);

        if (otherAttributes != null && !otherAttributes.isEmpty()) {
            for (final SAMRecord.SAMTagAndValue attribute : otherAttributes) {
                Utils.nonNull(attribute, "other attributes contain null attributes");
                result.setAttribute(attribute.tag, attribute.value);
            }
        }

        if (mismatches != NO_NM) {
            result.setAttribute(SAMTag.NM.name(), mismatches);
        }
        if (alnScore != NO_AS) {
            result.setAttribute(SAMTag.AS.name(), alnScore);
        }
        return result;
    }

    /**
     * Returns a cigar where all clips are transformed into either soft or hard clips.
     *
     * <p>
     *     After such transformation, any resulting adjacent clips are folded into a single element.
     * </p>
     * @param cigar the input cigar.
     * @param clipOperator the output cigar unified clipping operator (must be either {@link CigarOperator#S} or {@link CigarOperator#H}.
     *
     * @throws IllegalArgumentException if the input cigar cannot be {@code null}, or if it is invalid.
     * @return the output cigar after the transformation.
     */
    @VisibleForTesting
    static Cigar softOrHardReclip(final Cigar cigar, final CigarOperator clipOperator) {
        Utils.nonNull(cigar, "the input cigar cannot be null");
        final List<CigarElement> elements = cigar.getCigarElements();
        final int elementsSize = elements.size();
        if (elementsSize < 1) { // the only operand cannot be a clip.
            return new Cigar();
        } else if (elements.get(0).getOperator().isClipping() || elements.get(elementsSize - 1).getOperator().isClipping()) {
            final List<CigarElement> resultElements = new ArrayList<>(elements.size());
            // Stack-like construction of the result-elements list
            // where we fold-left the clips elements that we encounter
            CigarElement lastElement = null; // caches the last cigar element added.
            for (final CigarElement element : elements) {
                final CigarOperator operator = element.getOperator();
                if (operator.isClipping()) {
                    if (lastElement != null && lastElement.getOperator().isClipping()) { // fold-left merging clipping operations.
                        final int newLength = element.getLength() + lastElement.getLength();
                        resultElements.set(resultElements.size() - 1, lastElement = new CigarElement(newLength, clipOperator));
                    } else if (operator != clipOperator) { // cannot reuse the element instance as it has the "wrong" operator.
                        resultElements.add(lastElement = new CigarElement(element.getLength(), clipOperator));
                    } else { // has the right clip-operator and is the first in the run of clipping operations if any.
                        resultElements.add(lastElement = element);
                    }
                } else { // is not a clipping, just add it as it is.
                    resultElements.add(lastElement = element);
                }
            }
            return new Cigar(resultElements);
        } else {
            return new Cigar(elements);
        }
    }
}
