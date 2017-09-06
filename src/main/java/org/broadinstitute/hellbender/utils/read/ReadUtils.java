package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import htsjdk.samtools.cram.build.CramIO;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.recalibration.EventType;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;

/**
 * A miscellaneous collection of utilities for working with reads, headers, etc.
 * Static methods only, please.
 */
public final class ReadUtils {

    private ReadUtils() {
    }

    private static final Logger logger = LogManager.getLogger();

    /**
     * The default quality score for an insertion or deletion, if
     * none are provided for this read.
     */
    public static final byte DEFAULT_INSERTION_DELETION_QUAL = (byte)45;

    // Base Quality Score Recalibrator specific attribute tags
    public static final String BQSR_BASE_INSERTION_QUALITIES = "BI";                // base qualities for insertions
    public static final String BQSR_BASE_DELETION_QUALITIES = "BD";                 // base qualities for deletions

    public static final int CLIPPING_GOAL_NOT_REACHED = -1;
    private static final int DEFAULT_ADAPTOR_SIZE = 100;

    public static final String ORIGINAL_BASE_QUALITIES_TAG = SAMTag.OQ.name();

    /**
     * BAM file magic value that starts every bam file
     */
    public static final byte[] BAM_MAGIC = "BAM\1".getBytes();

    /**
     * HACK: This is used to make a copy of a read.
     * Really, SAMRecord should provide a copy constructor or a factory method.
     */
    public static SAMRecord cloneSAMRecord(final SAMRecord originalRead) {
        if (originalRead == null) {
            return null;
        }
        try {
            return (SAMRecord)originalRead.clone();
        } catch (final CloneNotSupportedException e) {
            throw new IllegalStateException(e);
        }
    }

    /**
     * HACK: This is used to make a copy of a header.
     * Really, SAMFileHeader should provide a copy constructor or a factory method.
     */

    public static SAMFileHeader cloneSAMFileHeader( final SAMFileHeader header ) {
        if (header == null) return null;
        return header.clone();
    }

    /**
     * Checks whether read is a headerless SAMRecordToGATKReadAdapter, and if it is, sets its
     * header to the provided header.
     *
     * @param read A potentially headerless GATKRead
     * @param header header to store in the read, if it's a headerless SAMRecord-backed read
     */
    public static void restoreHeaderIfNecessary( final GATKRead read, final SAMFileHeader header ) {
        if ( read instanceof SAMRecordToGATKReadAdapter ) {
            SAMRecordToGATKReadAdapter readAdapter = (SAMRecordToGATKReadAdapter)read;
            if ( ! readAdapter.hasHeader() ) {
                readAdapter.setHeader(header);
            }
        }
    }

    /**
     * Retrieve the original base qualities of the given read, if present,
     * as stored in the OQ attribute.
     *
     * @param read read to check
     * @return original base qualities as stored in the OQ attribute, or null
     *         if the OQ attribute is not present
     */
    public static byte[] getOriginalBaseQualities( final GATKRead read ) {
        if ( ! read.hasAttribute(ORIGINAL_BASE_QUALITIES_TAG) ) {
            return null;
        }
        final String oqString = read.getAttributeAsString(ORIGINAL_BASE_QUALITIES_TAG);
        return !oqString.isEmpty() ? SAMUtils.fastqToPhred(oqString) : null;
    }

    /**
     * Returns the base qualities for the read as a string.
     *
     * @param read read whose base qualities should be returned
     * @return Base qualities string as printable ASCII values (encoded as a FASTQ string).
     */
    public static String getBaseQualityString( final GATKRead read ) {
        Utils.nonNull(read);
        if ( Arrays.equals(SAMRecord.NULL_QUALS, read.getBaseQualities()) ) {
            return SAMRecord.NULL_QUALS_STRING;
        }
        return SAMUtils.phredToFastq(read.getBaseQualities());
    }

    /**
     * Set the base qualities from a string of ASCII encoded values
     * @param read read whose base qualities should be set
     * @param baseQualityString ASCII encoded (encoded as a FASTQ string) values of base qualities.
     */
    public static void setBaseQualityString(final GATKRead read, final String baseQualityString) {
        Utils.nonNull(read);
        Utils.nonNull(baseQualityString);

        if ( SAMRecord.NULL_QUALS_STRING.equals(baseQualityString) ) {
            read.setBaseQualities(SAMRecord.NULL_QUALS);
        } else {
            read.setBaseQualities(SAMUtils.fastqToPhred(baseQualityString));
        }
    }

    /**
     * Returns the reference index in the given header of the read's contig,
     * or {@link SAMRecord#NO_ALIGNMENT_REFERENCE_INDEX} if the read is unmapped.
     *
     * @param read read whose reference index to look up
     * @param header SAM header defining contig indices
     * @return the reference index in the given header of the read's contig,
     *         or {@link SAMRecord#NO_ALIGNMENT_REFERENCE_INDEX} if the read is unmapped.
     */
    public static int getReferenceIndex( final GATKRead read, final SAMFileHeader header ) {
        if ( read.isUnmapped() ) {
            return SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
        }

        return header.getSequenceIndex(read.getContig());
    }

    /**
     * Returns the reference index in the given header of the read's assigned contig.
     *
     * Unlike {@link #getReferenceIndex}, which returns {@link SAMRecord#NO_ALIGNMENT_REFERENCE_INDEX}
     * for all unmapped reads, this method will return a reference index for unmapped
     * reads that are assigned a nominal position (eg., the position of their mates),
     * which is useful for sorting.
     *
     * @param read read whose reference index to look up
     * @param header SAM header defining contig indices
     * @return the reference index in the given header of the read's contig,
     *         or {@link SAMRecord#NO_ALIGNMENT_REFERENCE_INDEX} if the read's contig
     *         is not found in the header
     */
    public static int getAssignedReferenceIndex( final GATKRead read, final SAMFileHeader header ) {
        return header.getSequenceIndex(read.getAssignedContig());
    }

    /**
     * Checks whether the provided read has an assigned position. This is different than checking
     * unmapped status, since unmapped reads are often assigned a nominal position (eg., the position
     * of their mapped mate). A read is considered to have no assigned position if its assigned contig
     * is either {@code null} or {@link ReadConstants#UNSET_CONTIG}, or its assigned start position is
     * {@link ReadConstants#UNSET_POSITION}, regardless of whether the read is actually marked as mapped
     * or unmapped.
     *
     * @param read read to check
     * @return true if the read has no assigned position, otherwise false
     */
    public static boolean readHasNoAssignedPosition( final GATKRead read ) {
        // Check actual assigned positions rather than unmapped status, so that unmapped reads with
        // assigned positions will be considered to have a position
        return read.getAssignedContig() == null ||
               read.getAssignedContig().equals(ReadConstants.UNSET_CONTIG) ||
               read.getAssignedStart() == ReadConstants.UNSET_POSITION;
    }

    /**
     * Returns the reference index in the given header of the contig of the read's mate,
     * or {@link SAMRecord#NO_ALIGNMENT_REFERENCE_INDEX} if the read's mate is unmapped.
     *
     * @param read read whose mate's reference index to look up
     * @param header SAM header defining contig indices
     * @return the reference index in the given header of the contig of the read's mate,
     *         or {@link SAMRecord#NO_ALIGNMENT_REFERENCE_INDEX} if the read's mate is unmapped.
     */
    public static int getMateReferenceIndex( final GATKRead read, final SAMFileHeader header ) {
        if ( read.mateIsUnmapped() ) {
            return SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
        }

        return header.getSequenceIndex(read.getMateContig());
    }

    /**
     * Returns a {@link SAMReadGroupRecord} object corresponding to the provided read's read group.
     *
     * @param read read whose read group to retrieve
     * @param header SAM header containing read groups
     * @return a {@link SAMReadGroupRecord} object corresponding to the provided read's read group,
     *         or null if the read has no read group
     */
    public static SAMReadGroupRecord getSAMReadGroupRecord( final GATKRead read, final SAMFileHeader header ) {
        final String readGroupName = read.getReadGroup();
        return readGroupName != null ? header.getReadGroup(readGroupName) : null;
    }

    /**
     * Returns the platform associated with the provided read's read group.
     *
     * @param read read whose platform information to retrieve
     * @param header SAM header containing read groups
     * @return the platform for the provided read's read group as a String,
     *         or null if the read has no read group.
     */
    public static String getPlatform( final GATKRead read, final SAMFileHeader header ) {
        final SAMReadGroupRecord readGroup = getSAMReadGroupRecord(read, header);
        return readGroup != null ? readGroup.getPlatform() : null;
    }

    /**
     * Returns the platform unit associated with the provided read's read group.
     *
     * @param read read whose platform unit to retrieve
     * @param header SAM header containing read groups
     * @return the platform unit for the provided read's read group as a String,
     *         or null if the read has no read group.
     */
    public static String getPlatformUnit( final GATKRead read, final SAMFileHeader header ) {
        final SAMReadGroupRecord readGroup = getSAMReadGroupRecord(read, header);
        return readGroup != null ? readGroup.getPlatformUnit() : null;
    }

    /**
     * Returns the library associated with the provided read's read group.
     *
     * @param read read whose library to retrieve
     * @param header SAM header containing read groups
     * @return the library for the provided read's read group as a String,
     *         or null if the read has no read group.
     */
    public static String getLibrary( final GATKRead read, final SAMFileHeader header ) {
        final SAMReadGroupRecord readGroup = getSAMReadGroupRecord(read, header);
        return readGroup != null ? readGroup.getLibrary() : null;
    }

    /**
     * Returns the sample name associated with the provided read's read group.
     *
     * @param read read whose sample name to retrieve
     * @param header SAM header containing read groups
     * @return the sample name for the provided read's read group as a String,
     *         or null if the read has no read group.
     */
    public static String getSampleName( final GATKRead read, final SAMFileHeader header ) {
        final SAMReadGroupRecord readGroup = getSAMReadGroupRecord(read, header);
        return readGroup != null ? readGroup.getSample() : null;
    }

    /**
     * Returns the read's unclipped start if the read is on the forward strand,
     * or the read's unclipped end if the read is on the reverse strand.
     *
     * @param read read whose stranded unclipped start to retrieve
     * @return the read's unclipped start if the read is on the forward strand,
     *         or the read's unclipped end if the read is on the reverse strand.
     */
    public static int getStrandedUnclippedStart( final GATKRead read ) {
        return read.isReverseStrand() ? read.getUnclippedEnd() : read.getUnclippedStart();
    }

    public static boolean isEmpty(final SAMRecord read) {
        return read.getReadBases() == null || read.getReadLength() == 0;
    }

    public static String prettyPrintSequenceRecords( final SAMSequenceDictionary sequenceDictionary ) {
        final String[] sequenceRecordNames = new String[sequenceDictionary.size()];
        int sequenceRecordIndex = 0;
        for (final SAMSequenceRecord sequenceRecord : sequenceDictionary.getSequences()) {
            sequenceRecordNames[sequenceRecordIndex++] = sequenceRecord.getSequenceName();
        }
        return Arrays.deepToString(sequenceRecordNames);
    }

    /**
     * @param read read to check
     * @return true if the read is paired and has a mapped mate, otherwise false
     */
    public static boolean readHasMappedMate( final GATKRead read ) {
        return read.isPaired() && ! read.mateIsUnmapped();
    }

    /**
     * Check whether the given String represents a legal attribute name according to the SAM spec,
     * and throw an exception if it doesn't.
     *
     * Legal attribute names are two characters long, start with a letter, and end with a letter or digit.
     *
     * @param attributeName name to check
     * @throws IllegalArgumentException if the attribute name is illegal according to the SAM spec.
     */
    public static void assertAttributeNameIsLegal( final String attributeName ) {
        if ( attributeName == null ||
             attributeName.length() != 2 ||
             ! Character.isLetter(attributeName.charAt(0)) ||
             ! Character.isLetterOrDigit(attributeName.charAt(1)) ) {

            throw new IllegalArgumentException("Read attribute " + attributeName + " invalid: attribute names must be non-null two-character Strings matching the pattern /[A-Za-z][A-Za-z0-9]/");
        }
    }

    /**
     * Encapsulates a integer attribute into an {@link OptionalInt} instance.
     * @param read the input read.
     * @param tag the attribute tag name.
     * @throws IllegalArgumentException if {@code read} or {@code tag} are {@code null}.
     * @throws GATKException.ReadAttributeTypeMismatch if the value provided for that attribute is not an integer.
     * @return never {@code null}, but perhaps empty indicating that no value was provided for this attribute.
     */
    public static OptionalInt getOptionalIntAttribute(final SAMRecord read, final String tag) {
        Utils.nonNull(read);
        Utils.nonNull(tag);
        final Object obj = read.getAttribute(tag);
        if (obj == null) {
            return OptionalInt.empty();
        } else if (obj instanceof Number) {
            final Number num = (Number) obj;
            if (num.intValue() == num.doubleValue()) {
                return OptionalInt.of(num.intValue());
            } else {
                throw new GATKException.ReadAttributeTypeMismatch(String.format("no an integer value for record %s @ %s and tag %s: %s",
                        read.getReadName(), locationString(read), tag, obj));
            }
        } else {
            final String str = "" + obj;
            try {
                return OptionalInt.of(Integer.parseInt(str));
            } catch (final NumberFormatException ex) {
                throw new GATKException.ReadAttributeTypeMismatch(String.format("no an integer value for record %s @ %s and tag %s: %s",
                        read.getReadName(), locationString(read), tag, obj), ex);
            }
        }
    }

    /**
     * Encapsulates a integer attribute into an {@link OptionalInt} instance.
     * @param read the input read.
     * @param tag the attribute tag name.
     * @throws IllegalArgumentException if {@code read} or {@code tag} are {@code null}.
     * @throws GATKException.ReadAttributeTypeMismatch if the value provided for that attribute is not an integer.
     * @return never {@code null}, but perhaps empty indicating that no value was provided for this attribute.
     */
    public static OptionalInt getOptionalIntAttribute(final GATKRead read, final String tag) {
        Utils.nonNull(read);
        Utils.nonNull(tag);
        final Integer obj = read.getAttributeAsInteger(tag);
        return obj == null ? OptionalInt.empty() : OptionalInt.of(obj);
    }

    private static String locationString(final SAMRecord read) {
        final String refName = read.getReferenceName();
        if (refName == null || SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(refName)) {
            return SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
        } else {
            final int start = read.getAlignmentStart();
            return refName + ":" + start;
        }
    }

    /**
     * A marker to tell which end of the read has been clipped
     */
    public enum ClippingTail {
        LEFT_TAIL,
        RIGHT_TAIL
    }

    public static final int SAM_READ_PAIRED_FLAG = 0x1;
    public static final int SAM_PROPER_PAIR_FLAG = 0x2;
    public static final int SAM_READ_UNMAPPED_FLAG = 0x4;
    public static final int SAM_MATE_UNMAPPED_FLAG = 0x8;
    public static final int SAM_READ_STRAND_FLAG = 0x10;
    public static final int SAM_MATE_STRAND_FLAG = 0x20;
    public static final int SAM_FIRST_OF_PAIR_FLAG = 0x40;
    public static final int SAM_SECOND_OF_PAIR_FLAG = 0x80;
    public static final int SAM_NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    public static final int SAM_READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    public static final int SAM_DUPLICATE_READ_FLAG = 0x400;
    public static final int SAM_SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;

    /**
     * Construct a set of SAM bitwise flags from a GATKRead
     *
     * @param read read from which to construct the flags
     * @return SAM-compliant set of bitwise flags reflecting the properties in the given read
     */
    public static int getSAMFlagsForRead( final GATKRead read ) {
        int samFlags = 0;

        if ( read.isPaired() ) {
            samFlags |= SAM_READ_PAIRED_FLAG;
        }
        if ( read.isProperlyPaired() ) {
            samFlags |= SAM_PROPER_PAIR_FLAG;
        }
        if ( read.isUnmapped() ) {
            samFlags |= SAM_READ_UNMAPPED_FLAG;
        }
        if ( read.isPaired() && read.mateIsUnmapped() ) {
            samFlags |= SAM_MATE_UNMAPPED_FLAG;
        }
        if ( !read.isUnmapped() && read.isReverseStrand() ) {
            samFlags |= SAM_READ_STRAND_FLAG;
        }
        if ( read.isPaired() && ! read.mateIsUnmapped() && read.mateIsReverseStrand() ) {
            samFlags |= SAM_MATE_STRAND_FLAG;
        }
        if ( read.isFirstOfPair() ) {
            samFlags |= SAM_FIRST_OF_PAIR_FLAG;
        }
        if ( read.isSecondOfPair() ) {
            samFlags |= SAM_SECOND_OF_PAIR_FLAG;
        }
        if ( read.isSecondaryAlignment() ) {
            samFlags |= SAM_NOT_PRIMARY_ALIGNMENT_FLAG;
        }
        if ( read.failsVendorQualityCheck() ) {
            samFlags |= SAM_READ_FAILS_VENDOR_QUALITY_CHECK_FLAG;
        }
        if ( read.isDuplicate() ) {
            samFlags |= SAM_DUPLICATE_READ_FLAG;
        }
        if ( read.isSupplementaryAlignment() ) {
            samFlags |= SAM_SUPPLEMENTARY_ALIGNMENT_FLAG;
        }

        return samFlags;
    }

    /**
     * Finds the adaptor boundary around the read and returns the first base inside the adaptor that is closest to
     * the read boundary. If the read is in the positive strand, this is the first base after the end of the
     * fragment (Picard calls it 'insert'), if the read is in the negative strand, this is the first base before the
     * beginning of the fragment.
     *
     * There are two cases we need to treat here:
     *
     * 1) Our read is in the reverse strand :
     *
     *     <----------------------| *
     *   |--------------------->
     *
     *   in these cases, the adaptor boundary is at the mate start (minus one)
     *
     * 2) Our read is in the forward strand :
     *
     *   |---------------------->   *
     *     <----------------------|
     *
     *   in these cases the adaptor boundary is at the start of the read plus the inferred insert size (plus one)
     *
     * @param read the read being tested for the adaptor boundary
     * @return the reference coordinate for the adaptor boundary (effectively the first base IN the adaptor, closest to the read.
     * CANNOT_COMPUTE_ADAPTOR_BOUNDARY if the read is unmapped or the mate is mapped to another contig.
     */
    public static int getAdaptorBoundary(final GATKRead read) {
        if ( ! hasWellDefinedFragmentSize(read) ) {
            return CANNOT_COMPUTE_ADAPTOR_BOUNDARY;
        } else if ( read.isReverseStrand() ) {
            return read.getMateStart() - 1;           // case 1 (see header)
        } else {
            final int insertSize = Math.abs(read.getFragmentLength());    // the inferred insert size can be negative if the mate is mapped before the read (so we take the absolute value)
            return read.getStart() + insertSize;  // case 2 (see header)
        }
    }

    public static int CANNOT_COMPUTE_ADAPTOR_BOUNDARY = Integer.MIN_VALUE;

    /**
     * Can the adaptor sequence of read be reliably removed from the read based on the alignment of
     * read and its mate?
     *
     * @param read the read to check
     * @return true if it can, false otherwise
     */
    public static boolean hasWellDefinedFragmentSize(final GATKRead read) {
        if ( read.getFragmentLength() == 0 )
            // no adaptors in reads with mates in another chromosome or unmapped pairs
        {
            return false;
	    }
        if ( ! read.isPaired() )
            // only reads that are paired can be adaptor trimmed
        {
            return false;
	}
        if ( read.isUnmapped() || read.mateIsUnmapped() )
            // only reads when both reads are mapped can be trimmed
        {
            return false;
	}
//        if ( ! read.isProperlyPaired() )
//            // note this flag isn't always set properly in BAMs, can will stop us from eliminating some proper pairs
//            // reads that aren't part of a proper pair (i.e., have strange alignments) can't be trimmed
//            return false;
        if ( read.isReverseStrand() == read.mateIsReverseStrand() )
            // sanity check on isProperlyPaired to ensure that read1 and read2 aren't on the same strand
	    {
            return false;
        }

        if ( read.isReverseStrand() ) {
            // we're on the negative strand, so our read runs right to left
            return read.getEnd() > read.getMateStart();
        } else {
            // we're on the positive strand, so our mate should be to our right (his start + insert size should be past our start)
            return read.getStart() <= read.getMateStart() + read.getFragmentLength();
        }
    }

    /**
     * If a read starts in INSERTION, returns the first element length.
     *
     * Warning: If the read has Hard or Soft clips before the insertion this function will return 0.
     *
     * @param read
     * @return the length of the first insertion, or 0 if there is none (see warning).
     */
    public static int getFirstInsertionOffset(final GATKRead read) {
        final CigarElement e = read.getCigarElement(0);
        if ( e.getOperator() == CigarOperator.I ) {
            return e.getLength();
        } else {
            return 0;
        }
    }

    /**
     * If a read ends in INSERTION, returns the last element length.
     *
     * Warning: If the read has Hard or Soft clips after the insertion this function will return 0.
     *
     * @param read
     * @return the length of the last insertion, or 0 if there is none (see warning).
     */
    public static int getLastInsertionOffset(final GATKRead read) {
        final List<CigarElement> cigarElements = read.getCigarElements();
        final CigarElement e = cigarElements.get(cigarElements.size() - 1);
        if ( e.getOperator() == CigarOperator.I ) {
            return e.getLength();
        } else {
            return 0;
        }
    }

    /**
     * Calculates the reference coordinate for the beginning of the read taking into account soft clips but not hard clips.
     *
     * Note: getUnclippedStart() adds soft and hard clips, this function only adds soft clips.
     *
     * @return the unclipped start of the read taking soft clips (but not hard clips) into account
     */
    public static int getSoftStart(final GATKRead read) {
        Utils.nonNull(read, "read");

        int softStart = read.getStart();
        for (final CigarElement cig : read.getCigarElements()) {
            final CigarOperator op = cig.getOperator();

            if (op == CigarOperator.SOFT_CLIP) {
                softStart -= cig.getLength();
            } else if (op != CigarOperator.HARD_CLIP) {
                break;
            }
        }
        return softStart;
    }

    /**
     * Calculates the reference coordinate for the end of the read taking into account soft clips but not hard clips.
     *
     * Note: getUnclippedEnd() adds soft and hard clips, this function only adds soft clips.
     *
     * @return the unclipped end of the read taking soft clips (but not hard clips) into account
     */
    public static int getSoftEnd(final GATKRead read) {
        Utils.nonNull(read, "read");

        boolean foundAlignedBase = false;
        int softEnd = read.getEnd();
        final List<CigarElement> cigs = read.getCigarElements();
        for (int i = cigs.size() - 1; i >= 0; --i) {
            final CigarElement cig = cigs.get(i);
            final CigarOperator op = cig.getOperator();

            if (op == CigarOperator.SOFT_CLIP){ // assumes the soft clip that we found is at the end of the aligned read
                softEnd += cig.getLength();
            } else if (op != CigarOperator.HARD_CLIP) {
                foundAlignedBase = true;
                break;
            }
        }
        if( !foundAlignedBase ) { // for example 64H14S, the soft end is actually the same as the alignment end
            softEnd = read.getEnd();
        }
        return softEnd;
    }

    public static int getReadCoordinateForReferenceCoordinateUpToEndOfRead(final GATKRead read, final int refCoord, final ClippingTail tail) {
        final int leftmostSafeVariantPosition = Math.max(getSoftStart(read), refCoord);
        return getReadCoordinateForReferenceCoordinate(getSoftStart(read), read.getCigar(), leftmostSafeVariantPosition, tail, false);
    }

    /**
     * Pre-processes the results of {@link #getReadCoordinateForReferenceCoordinate(int, Cigar, int, boolean)} to take care of
     * two corner cases:
     *
     * 1. If clipping the right tail (end of the read) getReadCoordinateForReferenceCoordinate and fall inside
     * a deletion return the base after the deletion. If clipping the left tail (beginning of the read) it
     * doesn't matter because it already returns the previous base by default.
     *
     * 2. If clipping the left tail (beginning of the read) getReadCoordinateForReferenceCoordinate and the
     * read starts with an insertion, and you're requesting the first read based coordinate, it will skip
     * the leading insertion (because it has the same reference coordinate as the following base).
     *
     * @return the read coordinate corresponding to the requested reference coordinate for clipping.
     */
    public static int getReadCoordinateForReferenceCoordinate(final GATKRead read, final int refCoord, final ClippingTail tail) {
        return getReadCoordinateForReferenceCoordinate(getSoftStart(read), read.getCigar(), refCoord, tail, false);
    }

    /**
     * Returns the read coordinate corresponding to the requested reference coordinate.
     *
     * WARNING: if the requested reference coordinate happens to fall inside or just before a deletion (or skipped region) in the read, this function
     * will return the last read base before the deletion (or skipped region). This function returns a
     * Pair(int readCoord, boolean fallsInsideOrJustBeforeDeletionOrSkippedRegion) so you can choose which readCoordinate to use when faced with
     * a deletion (or skipped region).
     *
     * SUGGESTION: Use getReadCoordinateForReferenceCoordinate(GATKSAMRecord, int, ClippingTail) instead to get a
     * pre-processed result according to normal clipping needs. Or you can use this function and tailor the
     * behavior to your needs.
     *
     * @param read
     * @param refCoord the requested reference coordinate
     * @return the read coordinate corresponding to the requested reference coordinate. (see warning!)
     */
    public static Pair<Integer, Boolean> getReadCoordinateForReferenceCoordinate(GATKRead read, int refCoord) {
        return getReadCoordinateForReferenceCoordinate(getSoftStart(read), read.getCigar(), refCoord, false);
    }

    public static int getReadCoordinateForReferenceCoordinate(final int alignmentStart, final Cigar cigar, final int refCoord, final ClippingTail tail, final boolean allowGoalNotReached) {
        final Pair<Integer, Boolean> result = getReadCoordinateForReferenceCoordinate(alignmentStart, cigar, refCoord, allowGoalNotReached);
        int readCoord = result.getLeft();

        // Corner case one: clipping the right tail and falls on deletion, move to the next
        // read coordinate. It is not a problem for the left tail because the default answer
        // from getReadCoordinateForReferenceCoordinate is to give the previous read coordinate.
        if (result.getRight() && tail == ClippingTail.RIGHT_TAIL) {
            readCoord++;
        }

        // clipping the left tail and first base is insertion, go to the next read coordinate
        // with the same reference coordinate. Advance to the next cigar element, or to the
        // end of the read if there is no next element.
        final CigarElement firstElementIsInsertion = readStartsWithInsertion(cigar);
        if (readCoord == 0 && tail == ClippingTail.LEFT_TAIL && firstElementIsInsertion != null) {
            readCoord = Math.min(firstElementIsInsertion.getLength(), cigar.getReadLength() - 1);
        }

        return readCoord;
    }

    private static Pair<Integer, Boolean> getReadCoordinateForReferenceCoordinate(final int alignmentStart, final Cigar cigar, final int refCoord, final boolean allowGoalNotReached) {
        int readBases = 0;
        int refBases = 0;
        boolean fallsInsideDeletionOrSkippedRegion = false;
        boolean endJustBeforeDeletionOrSkippedRegion = false;
        boolean fallsInsideOrJustBeforeDeletionOrSkippedRegion = false;

        final int goal = refCoord - alignmentStart;  // The goal is to move this many reference bases
        if (goal < 0) {
            if (allowGoalNotReached) {
                return new MutablePair<>(CLIPPING_GOAL_NOT_REACHED, false);
            } else {
                throw new GATKException("Somehow the requested coordinate is not covered by the read. Too many deletions?");
            }
        }
        boolean goalReached = refBases == goal;

        final Iterator<CigarElement> cigarElementIterator = cigar.getCigarElements().iterator();
        while (!goalReached && cigarElementIterator.hasNext()) {
            final CigarElement cigarElement = cigarElementIterator.next();
            int shift = 0;

            if (cigarElement.getOperator().consumesReferenceBases() || cigarElement.getOperator() == CigarOperator.SOFT_CLIP) {
                if (refBases + cigarElement.getLength() < goal) {
                    shift = cigarElement.getLength();
                } else {
                    shift = goal - refBases;
                }

                refBases += shift;
            }
            goalReached = refBases == goal;

            if (!goalReached && cigarElement.getOperator().consumesReadBases()) {
                readBases += cigarElement.getLength();
            }

            if (goalReached) {
                // Is this base's reference position within this cigar element? Or did we use it all?
                final boolean endsWithinCigar = shift < cigarElement.getLength();

                // If it isn't, we need to check the next one. There should *ALWAYS* be a next one
                // since we checked if the goal coordinate is within the read length, so this is just a sanity check.
                if (!endsWithinCigar && !cigarElementIterator.hasNext()) {
                    if (allowGoalNotReached) {
                        return new MutablePair<>(CLIPPING_GOAL_NOT_REACHED, false);
                    } else {
                        throw new GATKException(String.format("Reference coordinate corresponds to a non-existent base in the read. This should never happen -- check read with alignment start: %s  and cigar: %s", alignmentStart, cigar));
                    }
                }

                CigarElement nextCigarElement = null;

                // if we end inside the current cigar element, we just have to check if it is a deletion (or skipped region)
                if (endsWithinCigar) {
                    fallsInsideDeletionOrSkippedRegion = (cigarElement.getOperator() == CigarOperator.DELETION || cigarElement.getOperator() == CigarOperator.SKIPPED_REGION);
                }// if we end outside the current cigar element, we need to check if the next element is an insertion, deletion or skipped region.
                else {
                    nextCigarElement = cigarElementIterator.next();

                    // if it's an insertion, we need to clip the whole insertion before looking at the next element
                    if (nextCigarElement.getOperator() == CigarOperator.INSERTION) {
                        readBases += nextCigarElement.getLength();
                        if (!cigarElementIterator.hasNext()) {
                            if (allowGoalNotReached) {
                                return new MutablePair<>(CLIPPING_GOAL_NOT_REACHED, false);
                            } else {
                                throw new GATKException(String.format("Reference coordinate corresponds to a non-existent base in the read. This should never happen -- check read with alignment start: %s  and cigar: %s", alignmentStart, cigar));
                            }
                        }

                        nextCigarElement = cigarElementIterator.next();
                    }

                    // if it's a deletion (or skipped region), we will pass the information on to be handled downstream.
                    endJustBeforeDeletionOrSkippedRegion = (nextCigarElement.getOperator() == CigarOperator.DELETION || nextCigarElement.getOperator() == CigarOperator.SKIPPED_REGION);
                }

                fallsInsideOrJustBeforeDeletionOrSkippedRegion = endJustBeforeDeletionOrSkippedRegion || fallsInsideDeletionOrSkippedRegion;

                // If we reached our goal outside a deletion (or skipped region), add the shift
                if (!fallsInsideOrJustBeforeDeletionOrSkippedRegion && cigarElement.getOperator().consumesReadBases()) {
                    readBases += shift;
                }// If we reached our goal just before a deletion (or skipped region) we need
                    // to add the shift of the current cigar element but go back to it's last element to return the last
                    // base before the deletion (or skipped region) (see warning in function contracts)
                else if (endJustBeforeDeletionOrSkippedRegion && cigarElement.getOperator().consumesReadBases()) {
                    readBases += shift - 1;
                }// If we reached our goal inside a deletion (or skipped region), or just between a deletion and a skipped region,
                    // then we must backtrack to the last base before the deletion (or skipped region)
                else if (fallsInsideDeletionOrSkippedRegion ||
                        (endJustBeforeDeletionOrSkippedRegion && nextCigarElement.getOperator().equals(CigarOperator.N)) ||
                        (endJustBeforeDeletionOrSkippedRegion && nextCigarElement.getOperator().equals(CigarOperator.D))) {
                    readBases--;
                }
            }
        }

        if (!goalReached) {
            if (allowGoalNotReached) {
                return new MutablePair<>(CLIPPING_GOAL_NOT_REACHED, false);
            } else {
                throw new GATKException("Somehow the requested coordinate is not covered by the read. Alignment " + alignmentStart + " | " + cigar);
            }
        }

        return Pair.of(readBases, fallsInsideOrJustBeforeDeletionOrSkippedRegion);
    }

    /**
     * Is a base inside a read?
     *
     * @param read                the read to evaluate
     * @param referenceCoordinate the reference coordinate of the base to test
     * @return true if it is inside the read, false otherwise.
     */
    public static boolean isInsideRead(final GATKRead read, final int referenceCoordinate) {
        return referenceCoordinate >= read.getStart() && referenceCoordinate <= read.getEnd();
    }

    /**
     * @see #readStartsWithInsertion(Cigar, boolean) with ignoreClipOps set to true
     */
    public static CigarElement readStartsWithInsertion(final Cigar cigarForRead) {
        return readStartsWithInsertion(cigarForRead, true);
    }

    public static CigarElement readStartsWithInsertion(final Cigar cigarForRead, final boolean ignoreSoftClipOps) {
        return readStartsWithInsertion(cigarForRead.getCigarElements(), ignoreSoftClipOps);
    }

    /**
     * Checks if a read starts with an insertion.
     *
     * @param cigarElementsForRead    the CIGAR to evaluate
     * @param ignoreSoftClipOps   should we ignore S operators when evaluating whether an I operator is at the beginning?  Note that H operators are always ignored.
     * @return the element if it's a leading insertion or null otherwise
     */
    public static CigarElement readStartsWithInsertion(final List<CigarElement> cigarElementsForRead, final boolean ignoreSoftClipOps) {
        for ( final CigarElement cigarElement : cigarElementsForRead ) {
            if ( cigarElement.getOperator() == CigarOperator.INSERTION ) {
                return cigarElement;
            } else if ( cigarElement.getOperator() != CigarOperator.HARD_CLIP && ( !ignoreSoftClipOps || cigarElement.getOperator() != CigarOperator.SOFT_CLIP) ) {
                break;
            }
        }
        return null;
    }

    /**
     * Returns the reverse complement of the read bases
     *
     * @param bases the read bases
     * @return the reverse complement of the read bases
     */
    public static String getBasesReverseComplement(final byte[] bases) {
        String reverse = "";
        for (int i = bases.length-1; i >=0; i--) {
            reverse += (char) BaseUtils.getComplement(bases[i]);
        }
        return reverse;
    }

    /**
     * Returns the reverse complement of the read bases
     *
     * @param read the read
     * @return the reverse complement of the read bases
     */
    public static String getBasesReverseComplement(final GATKRead read) {
        return getBasesReverseComplement(read.getBases());
    }

    /**
     * Calculate the maximum read length from the given list of reads.
     * @param reads list of reads
     * @return      non-negative integer
     */
    public static int getMaxReadLength( final List<GATKRead> reads ) {
        if( reads == null ) { throw new IllegalArgumentException("Attempting to check a null list of reads."); }

        int maxReadLength = 0;
        for( final GATKRead read : reads ) {
            maxReadLength = Math.max(maxReadLength, read.getLength());
        }
        return maxReadLength;
    }

    /**
     * Creates an "empty" read with the provided read's read group and mate
     * information, but empty (not-null) fields:
     *  - Cigar String
     *  - Read Bases
     *  - Base Qualities
     *
     * Use this method if you want to create a new empty read based on
     * another read
     *
     * @param read a read to copy fields from
     * @return a read with no bases but safe for the GATK
     */
    public static GATKRead emptyRead( final GATKRead read ) {
        final GATKRead emptyRead = read.copy();

        emptyRead.setCigar("");
        emptyRead.setBases(new byte[0]);
        emptyRead.setBaseQualities(new byte[0]);

        emptyRead.clearAttributes();
        String readGroup = read.getReadGroup();
        if (readGroup != null) {
            emptyRead.setAttribute(SAMTag.RG.name(), readGroup);
        }

        return emptyRead;
    }

    public static void setInsertionBaseQualities( final GATKRead read, final byte[] quals) {
        read.setAttribute(BQSR_BASE_INSERTION_QUALITIES, quals == null ? null : SAMUtils.phredToFastq(quals));
    }

    public static void setDeletionBaseQualities( final GATKRead read, final byte[] quals) {
        read.setAttribute(BQSR_BASE_DELETION_QUALITIES, quals == null ? null : SAMUtils.phredToFastq(quals));
    }

    /**
     * @return whether or not this read has base insertion or deletion qualities (one of the two is sufficient to return true)
     */
    public static boolean hasBaseIndelQualities(final GATKRead read) {
        return read.hasAttribute(BQSR_BASE_INSERTION_QUALITIES) || read.hasAttribute(BQSR_BASE_DELETION_QUALITIES);
    }

    /**
     * @return the base deletion quality or null if read doesn't have one
     */
    public static byte[] getExistingBaseInsertionQualities(final GATKRead read) {
        return SAMUtils.fastqToPhred(read.getAttributeAsString(BQSR_BASE_INSERTION_QUALITIES));
    }

    /**
     * @return the base deletion quality or null if read doesn't have one
     */
    public static byte[] getExistingBaseDeletionQualities(final GATKRead read) {
        return SAMUtils.fastqToPhred( read.getAttributeAsString(BQSR_BASE_DELETION_QUALITIES));
    }

    /**
     * Default utility to query the base insertion quality of a read. If the read doesn't have one, it creates an array of default qualities (currently Q45)
     * and assigns it to the read.
     *
     * @return the base insertion quality array
     */
    public static byte[] getBaseInsertionQualities(final GATKRead read) {
        byte [] quals = getExistingBaseInsertionQualities(read);
        if( quals == null ) {
            quals = new byte[read.getBaseQualityCount()];
            Arrays.fill(quals, DEFAULT_INSERTION_DELETION_QUAL); // Some day in the future when base insertion and base deletion quals exist the samtools API will
            // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
        }
        return quals;
    }

    /**
     * Default utility to query the base deletion quality of a read. If the read doesn't have one, it creates an array of default qualities (currently Q45)
     * and assigns it to the read.
     *
     * @return the base deletion quality array
     */
    public static byte[] getBaseDeletionQualities(final GATKRead read) {
        byte[] quals = getExistingBaseDeletionQualities(read);
        if( quals == null ) {
            quals = new byte[read.getBaseQualityCount()];
            Arrays.fill(quals, DEFAULT_INSERTION_DELETION_QUAL);  // Some day in the future when base insertion and base deletion quals exist the samtools API will
            // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
        }
        return quals;
    }

    public static byte[] getBaseQualities( final GATKRead read, final EventType errorModel ) {
        switch( errorModel ) {
            case BASE_SUBSTITUTION:
                return read.getBaseQualities();
            case BASE_INSERTION:
                return getBaseInsertionQualities(read);
            case BASE_DELETION:
                return getBaseDeletionQualities(read);
            default:
                throw new GATKException("Unrecognized Base Recalibration type: " + errorModel );
        }
    }

    /**
     * Resets the quality scores of the reads to the orginal (pre-BQSR) ones.
     */
    public static GATKRead resetOriginalBaseQualities(final GATKRead read) {
        final byte[] originalQuals = ReadUtils.getOriginalBaseQualities(read);
        if ( originalQuals != null ){
            read.setBaseQualities(originalQuals);
        }
        return read;
    }

    /**
     * Check to ensure that the alignment makes sense based on the contents of the header.
     * @param header The SAM file header.
     * @param read The read to verify.
     * @return true if alignment agrees with header, false otherwise.
     */
    public static boolean alignmentAgreesWithHeader(final SAMFileHeader header, final GATKRead read) {
        final int referenceIndex = getReferenceIndex(read, header);
        // Read is aligned to nonexistent contig
        if( ! read.isUnmapped() && referenceIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX ) {
            return false;
        }
        final SAMSequenceRecord contigHeader = header.getSequence(referenceIndex);
        // Read is aligned to a point after the end of the contig
        return read.isUnmapped() || read.getStart() <= contigHeader.getSequenceLength();
    }

    /**
     * Create a common SAMFileWriter for use with GATK tools.
     *
     * @param outputFile - if this file has a .cram extension then a reference is required. Can not be null.
     * @param referenceFile - the reference source to use. Can not be null if a output file has a .cram extension.
     * @param header - header to be used for the output writer
     * @param preSorted - if true then the records must already be sorted to match the header sort order
     * @param createOutputBamIndex - if true an index will be created for .BAM and .CRAM files
     * @param createMD5 - if true an MD5 file will be created
     *
     * @return SAMFileWriter
     */
    public static SAMFileWriter createCommonSAMWriter(
            final File outputFile,
            final File referenceFile,
            final SAMFileHeader header,
            final boolean preSorted,
            boolean createOutputBamIndex,
            final boolean createMD5)
    {
        Utils.nonNull(outputFile);
        Utils.nonNull(header);

        if (createOutputBamIndex && header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            logger.warn("Skipping index file creation for: " +
                    outputFile.getAbsolutePath() +  ". Index file creation requires reads in coordinate sorted order.");
            createOutputBamIndex = false;
        }

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(createOutputBamIndex).setCreateMd5File(createMD5);
        return ReadUtils.createCommonSAMWriterFromFactory(factory, outputFile, referenceFile, header, preSorted);
    }

    /**
     * Create a common SAMFileWriter from a factory for use with GATK tools. Assumes that if the factory has been set
     * to create an index, the header must be set to coordinate sorted.
     *
     * @param outputFile if this file has a .cram extension then a reference is required. Can not be null.
     * @param referenceFile the reference source to use. Can not be null if a output file has a .cram extension.
     * @param header header to be used for the output writer
     * @param preSorted if true then records must already be sorted to match the header sort order
     * @param factory SAMFileWriterFactory factory to use
     * @return SAMFileWriter
     */
    public static SAMFileWriter createCommonSAMWriterFromFactory(
            final SAMFileWriterFactory factory,
            final File outputFile,
            final File referenceFile,
            final SAMFileHeader header,
            final boolean preSorted)
    {
        Utils.nonNull(outputFile);
        Utils.nonNull(header);

        if (null == referenceFile && outputFile.getName().endsWith(CramIO.CRAM_FILE_EXTENSION)) {
            throw new UserException("A reference file is required for writing CRAM files");
        }

        return factory.makeWriter(header.clone(), preSorted, outputFile, referenceFile);
    }

    /**
     * Validate that a file has CRAM contents by checking that it has a valid CRAM file header
     * (no matter what the extension).
     *
     * @param putativeCRAMFile File to check.
     * @return true if the file has a valid CRAM file header, otherwise false
     */
    public static boolean hasCRAMFileContents(final File putativeCRAMFile) {
        try (final FileInputStream fileStream = new FileInputStream(putativeCRAMFile);
             final BufferedInputStream bis = new BufferedInputStream(fileStream)) {
            return SamStreams.isCRAMFile(bis);
        }
        catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(e.getMessage());
        }
    }

    public static boolean isNonPrimary(GATKRead read) {
        return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped();
    }

    /**
     * is this base inside the adaptor of the read?
     *
     * There are two cases to treat here:
     *
     * 1) Read is in the negative strand => Adaptor boundary is on the left tail
     * 2) Read is in the positive strand => Adaptor boundary is on the right tail
     *
     * Note: We return false to all reads that are UNMAPPED or have an weird big insert size (probably due to mismapping or bigger event)
     *
     * @param read the read to test
     * @param basePos base position in REFERENCE coordinates (not read coordinates)
     * @return whether or not the base is in the adaptor
     */
    public static boolean isBaseInsideAdaptor(final GATKRead read, long basePos) {
        final int adaptorBoundary = ReadUtils.getAdaptorBoundary(read);
        if (adaptorBoundary == CANNOT_COMPUTE_ADAPTOR_BOUNDARY || read.getFragmentLength() > DEFAULT_ADAPTOR_SIZE)
            return false;

        return read.isReverseStrand() ? basePos <= adaptorBoundary : basePos >= adaptorBoundary;
    }

    /**
     * Pull out the sample names from a SAMFileHeader
     *
     * note that we use a TreeSet so that they are sorted
     *
     * @param header  the sam file header
     * @return list of strings representing the sample names
     */
    public static Set<String> getSamplesFromHeader( final SAMFileHeader header ) {
        // get all of the unique sample names
        final Set<String> samples = new TreeSet<>();
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();

        for ( SAMReadGroupRecord readGroup : readGroups ) {
            final String sample = readGroup.getSample();
            if ( sample != null ) {
                samples.add(sample);
            }
        }

        return samples;
    }

    /**
     * Validate that the expected input sort order is either "unsorted", or that it
     * matches the actualSortOrder. If validation fails a UserException is thrown
     * unless assumeSorted is true.
     *
     * @param actualSortOrder the actual sort order of the input
     * @param expectedSortOrder the sort order expected for this context
     * @param sourceName the name of the read source for inclusion in error messages
     * @param assumeSorted if true, no exception is thrown when the actualSortOrder
     *                     doesn't match the expectedSortOrder. An error messsage is
     *                     logged instead
     * @return boolean indicating if the validation passed
     * @throws UserException if the expectedSortOrder is anything other than "unsorted"
     * and the actualSortOrder doesn't match expectedSortOrder and assumeSorted is false
     */
    public static boolean validateExpectedSortOrder(
            final SAMFileHeader.SortOrder actualSortOrder,
            final SAMFileHeader.SortOrder expectedSortOrder,
            final boolean assumeSorted,
            final String sourceName)
    {
        boolean isValid = true;
        if (expectedSortOrder != SAMFileHeader.SortOrder.unsorted &&
                actualSortOrder != expectedSortOrder) {
            final String message = String.format("Input \"%s\" has sort order \"%s\" but \"%s\" is required.",
                    sourceName,
                    actualSortOrder.name(),
                    expectedSortOrder.name());
            isValid = false;
            if (assumeSorted) {
                logger.warn(message + " Assuming it's properly sorted anyway.");
            }
            else {
                throw new UserException(
                        message +
                                "If you believe the file to be sorted correctly, use " +
                                StandardArgumentDefinitions.ASSUME_SORTED_LONG_NAME +
                                "=true"
                );
            }
        }

        return isValid;
    }

    /**
     * Returns the offset (0-based index) of the first base in the read that is aligned against the reference.
     * <p>
     *     In most cases for mapped reads, this is typically equal to the sum of the size of soft-clipping at the
     *     beginning of the alignment.
     * </p>
     * <p>
     *     It returns {@code -1} when there is no base against the reference (e.g. the read is unmapped or is all in
     *     an insertion). This method will also return {@code -1} if the read does not contain a CIGAR.
     * </p>
     * <p>
     *     Notice that this index makes reference to the offset of that first base in the array returned by {@link #getBases()}, If you
     *     are after the first base in the original unclipped and not reverse-complemented read, you must use
     *     {@link #getFirstAlignedReadPosition} instead.
     * </p>
     *
     * @return a number between {@code -1} and the read length-1.
     */
    public static int getFirstAlignedBaseOffset(final GATKRead read) {
        Utils.nonNull(read, "the input read cannot be null");
        if (read.isUnmapped()) {
            return -1;
        } else {
            int result = 0;
            for (final CigarElement ce : read.getCigarElements()) {
                final int length = ce.getLength();
                if (length != 0) { // paranoia test.
                    final CigarOperator co = ce.getOperator();
                    if (co.isAlignment()) {
                        return result;
                    } else if (co.consumesReadBases()) {
                        result += length;
                    }
                }
            }
            return -1;
        }
    }

    /**
     * Returns the first position (1-based) on the original read in 5to3 prime orientation that is aligned
     * against a base on the reference.
     * <p>
     *     It will return {@code -1} to indicate that there is no base aligned against the reference (e.g.
     *     the read is unmapped or it totally included in and insertion and or clips. This is also the case if
     *     the read has an empty cigar for whatever reason.
     * </p>
     * <p>
     *     If what you are after is the offset of the first aligned base in the array returned by {@link #getBases()}
     *     then you should be looking into {@link #getFirstAlignedBaseOffset(GATKRead)}.
     * </p>
     * @return {@code -1} or a number between {@code 1} and the length of the original read. Never {@code 0}.
     */
    public static int getFirstAlignedReadPosition(final GATKRead read) {
        Utils.nonNull(read, "the input read cannot be null");
        if (read.isUnmapped()) {
            return -1;
        } else {
            final List<CigarElement> cigarElements = read.getCigarElements();
            if (cigarElements.size() < 1) {
                return -1;
            } else if (!read.isReverseStrand()) {
                final int offset = getFirstAlignedBaseOffset(read);
                return offset == -1 ? -1 : offset + 1 + CigarUtils.countLeftHardClippedBases(read.getCigar());
            } else { // isReverseStrand() == true
                final int offset = getAfterLastAlignmentBaseOffset(read);
                return offset == -1 ? -1 : read.getLength() - offset + 1 + CigarUtils.countRightHardClippedBases(read.getCigar());
            }
        }
    }

    /**
     * Returns the offset (0-based) after the last base in the record that is aligned against the reference.
     *
     * <p>
     *     It returns {@code -1} if there is not such a base (e.g. the record is unmapped or is totally enclosed in an
     *      insertion or clip). Also it will return {@code -1} if it has an empty cigar.
     * </p>
     * <p>
     *     Notice that this does not take in consideration whether the record is mapped against the forward and backward
     *     strand. If what you want to find out is the last base position on the original unclipped and non-reverve complemented read sequence that is mapped
     *     against the reference the use {@link #getLastAlignedReadPosition(GATKRead)} instead.
     * </p>
     *
     * @return a number between {@code -1} and the read's length inclusive.
     */
    public static int getAfterLastAlignmentBaseOffset(final GATKRead read) {
        if (read.isUnmapped()) {
            return -1;
        } else {
            int result = read.getLength();
            final List<CigarElement> cigarElements = read.getCigarElements();
            final ListIterator<CigarElement> cigarElementIterator = cigarElements.listIterator(cigarElements.size());
            while (cigarElementIterator.hasPrevious()) {
                final CigarElement ce = cigarElementIterator.previous();
                final int length = ce.getLength();
                if (length != 0) { // paranohic test.
                    final CigarOperator co = ce.getOperator();
                    if (co.isAlignment()) {
                        return result;
                    } else if (co.consumesReadBases()) {
                        result -= length;
                    }
                }
            }
            return -1;
        }
    }

    /**
     * Returns the last position (1-based) on the original unclipped nor reverse-complemented read that is mapped against the reference.
     * <p>
     *     It will return {@code -1} if no position in the read is aligned against the reference (e.g. unmapped).
     * </p>
     *
     * @return {@code -1} or any number between 1 and the length of the original read (including hard-clips). Never {@code 0}.
     */
    public static int getLastAlignedReadPosition(final GATKRead read) {
        if (read.isUnmapped()) {
            return -1;
        } else {
            final List<CigarElement> cigarElements = read.getCigarElements();
            if (cigarElements.isEmpty()) {
                return -1;
            } else if (!read.isReverseStrand()) {
                final int offset = getAfterLastAlignmentBaseOffset(read);
                return offset == -1 ? -1 : offset + CigarUtils.countLeftHardClippedBases(read.getCigar());
            } else { // isReverseStrand
                final int offset = getFirstAlignedBaseOffset(read);
                return offset == -1 ? -1 : read.getLength() - offset + CigarUtils.countRightHardClippedBases(read.getCigar());
            }
        }
    }
}
