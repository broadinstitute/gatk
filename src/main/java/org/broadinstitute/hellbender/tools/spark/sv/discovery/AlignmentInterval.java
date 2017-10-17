package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.hadoop.yarn.webapp.hamlet.Hamlet;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.AlnModType;
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
import org.broadinstitute.hellbender.utils.report.GATKReportColumnFormat;

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
     * when {@code alnModType} is not {@link AlnModType#NONE}, fields "mapQual", "mismatches", "alnScore" should be
     * viewed with care as they were either simply copied from the original alignment (e.g. "mapQual") or
     * set to some special value (e.g. "mismatches"), which is wrong strictly speaking.
     * (We didn't recompute them because that would require expensive SW re-alignment.)
     */
    public final AlnModType alnModType;


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
        this.alnModType = AlnModType.NONE;
    }

    public AlignmentInterval(final String contig, final int start, final boolean forwardStrand, final Cigar cigar,
                             final int mapQual, final int mismatches, final int alnScore) {
        Utils.nonNull(contig, "the input contig name cannot be negative");
        Utils.nonNull(cigar, "the input cigar cannot be negative");
        ParamUtils.isPositive(start, "the start position must be positive");
        Utils.validateArg(mapQual == SAMRecord.UNKNOWN_MAPPING_QUALITY || mapQual >= 0, "the mapping quality must be positive");
        Utils.validateArg(mismatches == AlignmentInterval.NO_NM || mismatches >= 0, "the number of mismatches must be positive");
        Utils.validateArg(alnScore == AlignmentInterval.NO_AS || alnScore >= 0, "the alignment score must be positive");
        this.cigarAlong5to3DirectionOfContig = forwardStrand ? cigar : CigarUtils.invertCigar(cigar);
        this.referenceSpan = new SimpleInterval(contig, start, start + cigar.getReferenceLength() - 1);
        this.startInAssembledContig = 1 + CigarUtils.countLeftClippedBases(this.cigarAlong5to3DirectionOfContig);
        this.endInAssembledContig = CigarUtils.countUnclippedReadBases(cigar) - CigarUtils.countRightClippedBases(cigar);
        this.alnModType = AlnModType.NONE;
        this.alnScore = alnScore;
        this.mapQual = mapQual;
        this.mismatches = mismatches;
        this.forwardStrand = forwardStrand;
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
        this.alnModType = AlnModType.NONE;
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
        this.alnModType = AlnModType.NONE;
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
        this.alnModType = AlnModType.NONE;
    }

    public AlignmentInterval(final SimpleInterval referenceSpan, final int startInAssembledContig, final int endInAssembledContig,
                             final Cigar cigarAlong5to3DirectionOfContig, final boolean forwardStrand,
                             final int mapQual, final int mismatches, final int alignerScore, final AlnModType modType) {
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
        alnModType = AlnModType.values()[input.readInt()];
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

    public PositionMatcher getPositionMatcher() {
        return new PositionMatcher();
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
        result = 31 * result + alnModType.hashCode();
        return result;
    }

    /**
     * Composes the Supplementary Alignment string that corresponds to the information contain in this interval with respect to
     * how it maps against the reference.
     * <p>
     * The format is the one described in {@link #AlignmentInterval(String)}. Notice that no ';' is appended at the end.
     * <p>
     *
     * @return never {@code null}.
     */
    public String toSumpplementaryAlignmentString() {
        final StringBuilder builder = new StringBuilder(100);
        builder.append(referenceSpan.getContig()).append(',')
               .append(referenceSpan.getStart()).append(',')
               .append(forwardStrand ? '+' : '-').append(',')
               .append(forwardStrand
                       ? cigarAlong5to3DirectionOfContig
                       : CigarUtils.invertCigar(cigarAlong5to3DirectionOfContig)).append(',')
               .append(mapQual);
        if (mismatches != NO_NM) {
            builder.append(',').append(mismatches);
        }
        if (alnScore != NO_AS) {
            builder.append(',').append(alnScore);
        } else {
            builder.setLength(builder.length());
        }
        return builder.toString();
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
     */
    public SAMRecord toSAMRecord(final SAMFileHeader header, final String name,
                                 final byte[] unclippedBases, final boolean hardClip, final int otherFlags,
                                 final Collection<? extends SAMRecord.SAMTagAndValue> otherAttributes) {
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


    /**
     * Convenient component to resolve the correspondence between matching based in the sequence vs the reference
     * within the alignment interval.
     */
    private class PositionMatcher {

        /**
         * Cached reference to the list of cigar-elements:
         */
        private final CigarElement[] cigarElements;

        /**
         *
         */
        private final AlignmentInterval interval;

        /**
         * Contains the position of the sequence (1-based) of the first base in the first cigar element that consumes
         * bases sequence bases starting from the i-th element (0-based.
         * <p>
         * So if the i-th cigar element is an alignment, insertion or clipping then the value in the i-th position is the
         * same as the position in the seq of the first base overlapped by that cigar, quite straight forward.
         * </p>
         * <p>In contrast when the i-th cigar element is something like a deletion, padding or N operator one need to look
         * for the next cigar element that consumes sequence bases</p>
         * <p>
         * For example if the 2nd element is a deletion (D) and the 3rd is
         * a alignment (M), then {@code cigarElementStartPositionInSeq[2]} is equal to
         * the first matching base in 3th element as deletion
         * do not consume sequence bases.
         * </p>
         * <p>
         *     the last position refers always contain the position after the last matched in the sequence.
         * </p>
         */
        private final int[] cigarElementStartPositionInSeq;

        /**
         * Similarly to {@link #cigarElementStartPositionInSeq} but this contains the corresponding first reference
         * base.
         * <p>
         *     If the alignment interval sits on the forward-strand, then values are ascendent whereas if the interval
         *     maps to the negative strand values are descenent.
         * </p>
         * <p>
         *     Clipping cigar elements before the first matching operator ha a value equal to the first reference bases matched by
         *     the interval if on the forward strand, the last matched based in on the reverse strand.
         * </p>
         * <p>
         *     Clipping cigar elements after the last matching operator ha a value just 1 over the last matched reference bases if
         *     on the forward strand, 1 under the first matched reference base if on the reverse strand.
         * </p>
         * <p>
         *     The last position always contains the position after (or before for the reverse strand) the last matched in the
         *     reference.
         * </p>
         */
        private final int[] cigarElementStartPositionInRef;

        private PositionMatcher() {
            interval = AlignmentInterval.this;
            cigarElements = cigarAlong5to3DirectionOfContig.getCigarElements().toArray(new CigarElement[cigarAlong5to3DirectionOfContig.numCigarElements()]);
            cigarElementStartPositionInRef = new int[cigarElements.length + 1];
            cigarElementStartPositionInSeq = new int[cigarElements.length + 1];
            int previousSeqPosition = cigarElementStartPositionInSeq[0] = 1;
            int previousRefPosition = cigarElementStartPositionInRef[0] = forwardStrand ? referenceSpan.getStart() : referenceSpan.getEnd();
            final int endOfReferenceAlignmentPosition = forwardStrand ? referenceSpan.getEnd() + 1 : referenceSpan.getStart() - 1;
            for (int i = 0; i < cigarElements.length; i++) {
                final CigarElement cigarElement  = cigarElements[i];
                final CigarOperator cigarOperator = cigarElement.getOperator();
                final int cigarLength = cigarElement.getLength();
                previousSeqPosition = cigarElementStartPositionInSeq[i + 1]
                        = cigarOperator.consumesReadBases() || cigarOperator.isClipping()
                             ? cigarLength + previousSeqPosition
                             : previousSeqPosition;
                previousRefPosition = cigarElementStartPositionInRef[i + 1]
                        = previousSeqPosition <= startInAssembledContig || previousRefPosition == endOfReferenceAlignmentPosition
                        ? previousRefPosition  // if before any alignment or after any alignment, the we just propagate the value forward.
                        : cigarOperator.consumesReferenceBases()
                            ? (forwardStrand ? cigarLength + previousRefPosition : previousRefPosition - cigarLength)
                            : previousRefPosition;
            }
        }

        /**
         * Return true iff the enclosing alignment intervals don't overlap in either the reference or the contig
         * or they do and the overlapping position in the contig match the same ones in the reference.
         * @param other the other interval position matcher.
         */
        public boolean isCompatibleWith(final PositionMatcher other) {
            Utils.nonNull(other);
            final boolean refOverlap = referenceSpan.overlaps(other.interval.referenceSpan);
            final boolean seqOverlap = interval.overlapOnContig(other.interval);
            if (refOverlap != seqOverlap) { // if ref doesnt overlap then seq must not to be compatible.
                return false;
            } else if (!refOverlap) { // neither overlap then they are also compatible.
                return true;
            } else {
                final int maxStart = Math.max(startInAssembledContig, other.interval.startInAssembledContig);
                final int minEnd = Math.min(endInAssembledContig, other.interval.endInAssembledContig);
                final int thisStartIndex = findEnclosingElementIndex(maxStart, 0);
                final int otherStartIndex = other.findEnclosingElementIndex(maxStart, 0);
                int nextPosition = maxStart;
                for (int i = thisStartIndex, j = otherStartIndex; nextPosition <= minEnd;) {
                    if (cigarElements[i].getOperator().isAlignment() != other.cigarElements[j].getOperator().isAlignment()) {
                        return false;
                    } else if (cigarElements[i].getOperator().isAlignment()) {
                        final int thisRefPosition = nextPosition - cigarElementStartPositionInSeq[i] + cigarElementStartPositionInRef[j];
                        final int otherRefPosition = nextPosition - other.cigarElementStartPositionInSeq[j] + other.cigarElementStartPositionInRef[j];
                        if (thisRefPosition != otherRefPosition) {
                            return false;
                        }
                    }
                    nextPosition = Math.min(cigarElementStartPositionInSeq[i + 1], cigarElementStartPositionInSeq[j + 1]);
                    i = findEnclosingElementIndex(nextPosition, i);
                    j = findEnclosingElementIndex(nextPosition, j);
                }
                return true;
            }
        }

        private int findEnclosingElementIndex(final int position, final int from) {
            for (int i = from + 1; i < cigarElementStartPositionInSeq.length; i++) {
                if (cigarElementStartPositionInSeq[i] > position) {
                    return i - 1;
                }
            }
            return -1;
        }
    }

    public boolean overlapOnContig(final AlignmentInterval other) {
        Utils.nonNull(other);
        return other.endInAssembledContig >= startInAssembledContig
                && endInAssembledContig >= other.startInAssembledContig;
    }

    public boolean compatibleWith(final AlignmentInterval other) {
        return getPositionMatcher().isCompatibleWith(other.getPositionMatcher());
    }

    /**
     * Returns the reciprocal alignment interval given the enclosing contig name and reference sequence length.
     * @param contigName
     * @param referenceLength
     * @return never {@code null}
     */
    public AlignmentInterval reciprocal(final String contigName, final int referenceLength) {
        final SimpleInterval newReferenceSpan = new SimpleInterval(contigName, startInAssembledContig, endInAssembledContig);
        final int newStartInAssembledContig = referenceSpan.getStart();
        final int newEndInAssembledContig = referenceSpan.getEnd();
        final Cigar cigar = CigarUtils.reciprocal(cigarAlongReference(), referenceSpan.getStart() - 1,
                referenceLength - referenceSpan.getEnd(), CigarOperator.H);
        return new AlignmentInterval(newReferenceSpan, newStartInAssembledContig, newEndInAssembledContig, cigar, forwardStrand, mapQual, mismatches, alnScore, alnModType);
    }

}
