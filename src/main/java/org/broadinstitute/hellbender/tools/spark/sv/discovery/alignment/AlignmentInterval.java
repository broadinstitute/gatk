package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import scala.Tuple2;

import java.util.*;

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

    /**
     * String used to separate several alignment-intervals in a
     * SA Tag like string.
     */
    public static final String SA_TAG_INTERVAL_SEPARATOR_STR = ";";

    /**
     * String used to separate fields in an alignment-interval string
     * representation.
     */
    public static final String SA_TAG_FIELD_SEPARATOR = ",";

    /**
     * Represents the lack of a value for a field in a string representation.
     */
    public static final String NO_VALUE_STR = ".";

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

    /**
     * when {@code alnModType} is not {@link ContigAlignmentsModifier.AlnModType#NONE}, fields "mapQual", "mismatches", "alnScore" should be
     * viewed with care as they were either simply copied from the original alignment (e.g. "mapQual") or
     * set to some special value (e.g. "mismatches"), which is wrong strictly speaking.
     * (We didn't recompute them because that would require expensive SW re-alignment.)
     */
    public final ContigAlignmentsModifier.AlnModType alnModType;


    /**
     * Compose an alignment interval instance from a SAM supplementary alignment formatted string.
     * <p>
     *     The input string format is:
     *     <pre>
     *         contig-name,start,strand,cigar,mq(,nm(,as)?)?
     *     </pre>
     *     where {@code nm} (number of mismatches) and {@code as} (alignment score) might be absent. The strand
     *     is symbolized as '+' for the forward-strand and '-' for the reverse-strand.
     *     <p>An example with different degrees of completeness:</p>
     *     <pre>
     *         chr10,1241241,+,10S1313M45I14M100H,30,5,66
     *         chr10,1241241,+,10S1313M45I14M100H,30,5
     *         chr10,1241241,+,10S1313M45I14M100H,30
     *     </pre>
     * </p>
     * <p>
     *     Additionally for mapping quality, number of mismatches and alignment score one
     *     can also represent that any value is unknown using {@value #NO_VALUE_STR}. This way
     *     one can provide a concrete value for alignment score while the number of mismatches is
     *     remains undefined.
     * </p>
     * <p>
     *     Example:
     *     <pre>
     *         chr10,1241241,+,10S1313M45I14M100H,30,.,66
     *     </pre>
     * </p>
     * <p>
     *     For mapping quality this would be equivalent to use the special constant {@link SAMRecord#UNKNOWN_MAPPING_QUALITY} ({@value SAMRecord#UNKNOWN_MAPPING_QUALITY}).
     * </p>
     * @param samSAtagString the input string.
     * @throws IllegalArgumentException if {@code str} is {@code null} or it does not look like a valid
     *                                  SA string.
     */
    public AlignmentInterval(final String samSAtagString) {
        Utils.nonNull(samSAtagString, "input str cannot be null");
        final int firstSeparatorIndex = samSAtagString.indexOf(SA_TAG_INTERVAL_SEPARATOR_STR);
        final String terminationTrimmedString;
        if (firstSeparatorIndex < 0) {
            terminationTrimmedString = samSAtagString;
        } else if (firstSeparatorIndex == samSAtagString.length() - 1) {
            terminationTrimmedString = samSAtagString.substring(0, firstSeparatorIndex);
        } else { // there are several intervals (separated by ;) in the input string.
            throw new IllegalArgumentException("the input string cannot contain more than one interval");
        }

        final String[] parts = terminationTrimmedString.split(SA_TAG_FIELD_SEPARATOR);
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
        final int mappingQuality = ParamUtils.inRange(parseZeroOrPositiveInt(parts[nextPartsIndex++], SAMRecord.UNKNOWN_MAPPING_QUALITY, "invalid mapping quality"), 0, 255, "the mapping quality must be in the range [0, 255]");
        final int mismatches = parts.length > nextPartsIndex ? parseZeroOrPositiveInt(parts[nextPartsIndex++], NO_NM, "invalid number of mismatches") : NO_NM;
        final int alignmentScore = parts.length > nextPartsIndex ? parseZeroOrPositiveInt(parts[nextPartsIndex], NO_AS, "invalid alignment score") : NO_AS;

        this.referenceSpan = new SimpleInterval(referenceContig, start,
                Math.max(start, cigar.getReferenceLength() + start - 1));
        this.startInAssembledContig = 1 + CigarUtils.countLeftClippedBases(cigar);
        this.endInAssembledContig = CigarUtils.countUnclippedReadBases(cigar) - CigarUtils.countRightClippedBases(cigar);
        this.mapQual = mappingQuality;
        this.mismatches = mismatches;
        this.alnScore = alignmentScore;
        this.forwardStrand = forwardStrand;
        this.cigarAlong5to3DirectionOfContig = cigar;
        this.alnModType = ContigAlignmentsModifier.AlnModType.NONE;
    }

    /**
     * Returns the SA tag string representation for this interval.
     * @return never {@code null}.
     */
    public String toSATagString() {
        return appendSATagString(new StringBuilder(100)).toString();
    }

    @Override
    public String toString() {
        return toSATagString();
    }

    /**
     * Appends the SA string representation of this interval into a builder.
     * @param builder where to append to.
     * @return reference to the input builder.
     * @throws IllegalArgumentException if the builder is {@code null}.
     */
    public StringBuilder appendSATagString(final StringBuilder builder) {
        Utils.nonNull(builder);
        builder.append(referenceSpan.getContig()).append(SA_TAG_FIELD_SEPARATOR)
               .append(referenceSpan.getStart()).append(SA_TAG_FIELD_SEPARATOR)
               .append(forwardStrand ? Strand.POSITIVE.toString() : Strand.NEGATIVE.toString()).append(SA_TAG_FIELD_SEPARATOR)
               .append(forwardStrand ? cigarAlong5to3DirectionOfContig.toString() : CigarUtils.invertCigar(cigarAlong5to3DirectionOfContig).toString()).append(SA_TAG_FIELD_SEPARATOR)
               .append(mapQual);
        if (mismatches != NO_NM || alnScore != NO_AS) {
            builder.append(SA_TAG_FIELD_SEPARATOR);
            if (mismatches == NO_NM) {
                builder.append(NO_VALUE_STR).append(SA_TAG_FIELD_SEPARATOR).append(alnScore);
            } else if (alnScore == NO_AS) {
                builder.append(mismatches);
            } else {
                builder.append(mismatches).append(SA_TAG_FIELD_SEPARATOR).append(alnScore);
            }
        }
        return builder;
    }

    private static int parseZeroOrPositiveInt(final String str, final int defaultValue, final String errorMessage) {
        if (str.equals(NO_VALUE_STR)) {
            return defaultValue;
        } else {
            try {
                return ParamUtils.isPositiveOrZero(Integer.parseInt(str), errorMessage + ": " + str);
            } catch (final NumberFormatException ex) {
                throw new IllegalArgumentException(errorMessage + "; not a valid integer in: " + str, ex);
            }
        }
    }

    /**
     * todo: we used to use htjdk Strand class, but it seems our custom Strand class in the SV package serves better here
     * Parses a strand text into a {@link Strand} object.
     *
     * @param strandText the input strand text.
     * @return never {@code null}, either {@link Strand#POSITIVE} or {@link Strand#NEGATIVE}.
     * @throws IllegalArgumentException if {@code strandText} is not a valid strand representation
     */
    private static Strand parseStrand(final String strandText) {
        try {
            return Strand.decode(strandText);
        } catch (final NoSuchElementException ex) {
            throw new IllegalArgumentException(ex);
        }
    }

    @VisibleForTesting
    public AlignmentInterval(final SAMRecord samRecord) {

        Utils.validateArg( ! Utils.nonNull(samRecord).getReadUnmappedFlag(),
                "sam record being used to construct AlignmentInterval is unmapped: " + samRecord.toString());

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
        this.alnModType = ContigAlignmentsModifier.AlnModType.NONE;
    }

    /**
     * Constructs an alignment interval that reflects on the mapping properties of a {@link GATKRead} instance.
     *
     * @param read the target read.
     * @throws IllegalArgumentException if {@code read} is {@code null}.
     */
    public AlignmentInterval(final GATKRead read) {

        Utils.validateArg( ! Utils.nonNull(read).isUnmapped(),
                "read being used to construct AlignmentInterval is unmapped: " + read.toString());

        final boolean isMappedReverse = read.isReverseStrand();
        this.referenceSpan = new SimpleInterval(read);
        this.startInAssembledContig = ReadUtils.getFirstAlignedReadPosition(read);
        this.endInAssembledContig = ReadUtils.getLastAlignedReadPosition(read);
        this.cigarAlong5to3DirectionOfContig = isMappedReverse ? CigarUtils.invertCigar(read.getCigar()) : read.getCigar();
        this.forwardStrand = !isMappedReverse;
        this.mapQual = read.getMappingQuality();
        this.mismatches = ReadUtils.getOptionalIntAttribute(read, SAMTag.NM.name()).orElse(NO_NM);
        this.alnScore = ReadUtils.getOptionalIntAttribute(read, SAMTag.AS.name()).orElse(NO_AS);
        this.alnModType = ContigAlignmentsModifier.AlnModType.NONE;
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
        this.alnModType = ContigAlignmentsModifier.AlnModType.NONE;
    }

    public AlignmentInterval(final SimpleInterval referenceSpan, final int startInAssembledContig, final int endInAssembledContig,
                             final Cigar cigarAlong5to3DirectionOfContig, final boolean forwardStrand,
                             final int mapQual, final int mismatches, final int alignerScore, final ContigAlignmentsModifier.AlnModType modType) {
        checkValidArgument(cigarAlong5to3DirectionOfContig, referenceSpan, startInAssembledContig, endInAssembledContig);

        this.referenceSpan = referenceSpan;
        this.startInAssembledContig = startInAssembledContig;
        this.endInAssembledContig = endInAssembledContig;

        this.cigarAlong5to3DirectionOfContig = cigarAlong5to3DirectionOfContig;

        this.forwardStrand = forwardStrand;
        this.mapQual = mapQual;
        this.mismatches = mismatches;
        this.alnScore = alignerScore;
        this.alnModType = modType;
    }

    @VisibleForTesting
    static final void checkValidArgument(final Cigar cigar, final SimpleInterval referenceSpan,
                                         final int readStart, final int readEnd) {

        final int softClippedBases = SvCigarUtils.checkCigarAndConvertTerminalInsertionToSoftClip(cigar).stream().filter(ce -> ce.getOperator().equals(CigarOperator.S)).mapToInt(CigarElement::getLength).sum();
        final int readLength = cigar.getReadLength() - softClippedBases;
        final int referenceLength = cigar.getReferenceLength();
        final boolean validState = referenceLength == referenceSpan.size() && readLength == (readEnd - readStart + 1);
        if ( ! validState) {
            throw new IllegalArgumentException("Encountering invalid arguments for constructing alignment,\t" +
                    "cigar: " + cigar.toString() + " ref.span: " + referenceSpan.toString() + " read span: " + readStart + "-" + readEnd);
        }
    }

    public boolean containsGapOfEqualOrLargerSize(final int gapSize) {
        return cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                .anyMatch(cigarElement ->
                        cigarElement.getOperator().isIndel() && cigarElement.getLength() >= gapSize);
    }

    public int getSizeOnRead() {
        return endInAssembledContig - startInAssembledContig + 1;
    }

    public boolean containsOnRef(final AlignmentInterval other) {
        return this.referenceSpan.contains(other.referenceSpan);
    }

    public boolean containsOnRead(final AlignmentInterval other) {
        return this.startInAssembledContig <= other.startInAssembledContig && this.endInAssembledContig >= other.endInAssembledContig;
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
    public static int overlapOnRefSpan(final AlignmentInterval one, final AlignmentInterval two) {

        if (!one.referenceSpan.getContig().equals(two.referenceSpan.getContig())) return 0;

        // dummy number for chr to be used in constructing SVInterval, since 2 input AI's both map to the same chr by this point
        final int dummyChr = 0;
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
        alnModType = ContigAlignmentsModifier.AlnModType.values()[input.readInt()];
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
        output.writeInt(alnModType.ordinal());
    }

    /**
     * Returns the cigar of this alignment interval along the reference.
     * <p>
     *     This is the same as {@link #cigarAlong5to3DirectionOfContig} if in
     *     the forward-strand, and its inverse if in the reverse-strand.
     * </p>
     * @return never {@code null}.
     */
    public Cigar cigarAlongReference() {
        return (forwardStrand) ? cigarAlong5to3DirectionOfContig : CigarUtils.invertCigar(cigarAlong5to3DirectionOfContig);
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
    public static final String PACKED_STRING_REP_SEPARATOR = "_";

    /**
     * @return A packed String representation of this alignment interval; intended for debugging or annotation usage
     * (both requires compactified message).
     */
    public String toPackedString() {
        return String.join(PACKED_STRING_REP_SEPARATOR, String.valueOf(startInAssembledContig),
                String.valueOf(endInAssembledContig), referenceSpan.toString(), (forwardStrand ? "+" : "-"),
                TextCigarCodec.encode(cigarAlong5to3DirectionOfContig),
                String.valueOf(mapQual), String.valueOf(mismatches), String.valueOf(alnScore),
                alnModType.toString());
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
        if (alnScore != that.alnScore) return false;
        if (!referenceSpan.equals(that.referenceSpan)) return false;
        if (!cigarAlong5to3DirectionOfContig.equals(that.cigarAlong5to3DirectionOfContig)) return false;
        return alnModType == that.alnModType;
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
        result = 31 * result + alnModType.ordinal();
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
        if (hardClip && SAMFlag.SECONDARY_ALIGNMENT.isUnset(otherFlags) && SAMFlag.SUPPLEMENTARY_ALIGNMENT.isUnset(otherFlags)) {
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

    /**
     * Given an alignment (represented by this AlignmentInterval)
     * and a reference span (in this case {@code otherRefSpan}),
     * the overlap between the two reference spans, if any,
     * has a corresponding span on the read to which the alignment belongs.
     * This method computes that 1-based, inclusive span on the read.
     *
     * The utility we have in mind, currently, is for
     * inferring intervals on a read
     * where it overlaps the requested {@code otherRefSpan}
     * so that one can extract corresponding the read sequence.
     *
     * If {@code otherRefSpan} does not overlap with the ref span of the alignment;
     * the computed value is a tuple 2 of (-1, -1).
     */
    public Tuple2<Integer, Integer> readIntervalAlignedToRefSpan(final SimpleInterval otherRefSpan) {
        final int start, end;
        if ( ! Utils.nonNull(otherRefSpan).overlaps(referenceSpan) ) { // ref spans doesn't overlap
            start = end = -1;
        } else if ( otherRefSpan.contains(referenceSpan) ) { // this ref span contained in other ref span, the entire read block by this alignment
            start = startInAssembledContig;
            end = endInAssembledContig;
        } else { // most complicated

            // computing because utility method requests hard clipped bases from the start of the read should NOT be counted
            final CigarElement firstCigarElement = cigarAlong5to3DirectionOfContig.getFirstCigarElement();
            final int hardClipOffset = firstCigarElement.getOperator().equals(CigarOperator.H) ? firstCigarElement.getLength() : 0;

            // computes how many bases on this reference span need to be walked, both from start and from end
            // examples:
            // this :   |------------|
            // that :           |---------------|
            // ovlap:           |----|
            // bases:   |<- 9 ->|
            // then walking distance from start would be 8=9-1, and from end would be 0
            // this :           |---------------|
            // that :   |------------|
            // ovlap:           |----|
            // bases:                |<-  12  ->|
            // then walking distance from start would be 0, and from end would be 11=12-1
            final int distOnRefForStart;
            final int distOnRefForEnd;
            if (forwardStrand) {
                distOnRefForStart = Math.max(0, otherRefSpan.getStart() - referenceSpan.getStart());
                distOnRefForEnd = Math.max(0, referenceSpan.getEnd() - otherRefSpan.getEnd());
            } else {
                distOnRefForStart = Math.max(0, referenceSpan.getEnd() - otherRefSpan.getEnd());
                distOnRefForEnd = Math.max(0, otherRefSpan.getStart() - referenceSpan.getStart());
            }
            // then using the walking distances on alignment's ref span to compute
            // associated walking distances on read interval (by definition must be a sub-interval of alignment's read span)
            int walkDistOnReadFromStart = 0;
            if ( distOnRefForStart != 0 ) {
                final int startPosOnRead = startInAssembledContig - hardClipOffset; // utility method requests no hard clipped counted
                walkDistOnReadFromStart = SvCigarUtils.computeAssociatedDistOnRead(cigarAlong5to3DirectionOfContig, startPosOnRead, distOnRefForStart, false);
            }
            int walkDistOnReadFromEnd = 0;
            if ( distOnRefForEnd != 0 ) {
                final int startPosOnRead = endInAssembledContig - hardClipOffset; // utility method requests no hard clipped counted
                walkDistOnReadFromEnd = SvCigarUtils.computeAssociatedDistOnRead(cigarAlong5to3DirectionOfContig, startPosOnRead, distOnRefForEnd, true);
            }

            start = startInAssembledContig + walkDistOnReadFromStart;
            end = endInAssembledContig - walkDistOnReadFromEnd;
        }
        return new Tuple2<>(start, end);
    }
}
