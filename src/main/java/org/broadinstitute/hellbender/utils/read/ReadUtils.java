package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import htsjdk.samtools.util.FileExtensions;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.recalibration.EventType;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.util.*;

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

    public static final int READ_INDEX_NOT_FOUND = -1;
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
     * @param read read to query
     * @return true if the read has a mate and that mate is mapped, otherwise false
     */
    public static boolean readHasMappedMate( final GATKRead read ) {
        return read.isPaired() && ! read.mateIsUnmapped();
    }

    /**
     * @param read read to query
     * @return true if the read has a mate and that mate is mapped, otherwise false
     */
    public static boolean readHasMappedMate( final SAMRecord read ) {
        return read.getReadPairedFlag() && ! read.getMateUnmappedFlag();
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
        } else if (obj instanceof Integer || obj instanceof Short) {
            final Number num = (Number) obj;
            return OptionalInt.of(num.intValue());
        } else if (obj instanceof CharSequence) {
            final String str = "" + obj;
            try {
                return OptionalInt.of(Integer.parseInt(str));
            } catch (final NumberFormatException ex) {
                throw new GATKException.ReadAttributeTypeMismatch(read, tag, "integer", ex);
            }
        } else {
            throw new GATKException.ReadAttributeTypeMismatch(read, tag, "integer", obj);
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

    /**
     * Helper method for interrogating if a read and its mate (if it exists) are unmapped
     * @param read a read with mate information to interrogate
     * @return true if this read and its are unmapped
     */
    public static boolean readAndMateAreUnmapped(GATKRead read) {
        return read.isUnmapped() && (!read.isPaired() || read.mateIsUnmapped());
    }

    /**
     * Interrogates the header to determine if the bam is expected to be sorted such that reads with the same name appear in order.
     * This can correspond to either a queryname sorted bam or a querygrouped bam (unordered readname groups)
     * @param header header corresponding to the bam file in question
     * @return true if the header has has the right readname group
     */
    public static boolean isReadNameGroupedBam(SAMFileHeader header) {
        return SAMFileHeader.SortOrder.queryname.equals(header.getSortOrder()) || SAMFileHeader.GroupOrder.query.equals(header.getGroupOrder());
    }

    /**
     * Create a map of reads overlapping {@code interval} to their mates by looking for all possible mates within some
     * maximum fragment size.  This is not guaranteed to find all mates, in particular near structural variant breakpoints
     * where mates may align far away.
     *
     * The algorithm is:
     * 1) make two maps of read name --> read for reads overlapping {@code interval}, one for first-of-pair reads and one
     *    for second-of-pair reads.
     * 2) For all reads in an expanded interval padded by {@code fragmentSize} on both sides look for a read of the same name
     *    that is second-of-pair if this read is first-of-pair or vice-versa.  If such a read is found then this is that read's mate.
     *
     * @param readsContext
     * @param fragmentSize the maximum distance on either side of {@code interval} to look for mates.
     * @return a map of reads ot their mates for all reads for which a mate could be found.
     */
    public static Map<GATKRead, GATKRead> getReadToMateMap(final ReadsContext readsContext, final int fragmentSize) {
        final Map<String, GATKRead> readOnes = new HashMap<>();
        final Map<String, GATKRead> readTwos = new HashMap<>();
        Utils.stream(readsContext.iterator()).forEach(read -> (read.isFirstOfPair() ? readOnes : readTwos).put(read.getName(), read));

        final Map<GATKRead, GATKRead> result = new HashMap<>();
        final SimpleInterval originalInterval = readsContext.getInterval();
        final SimpleInterval expandedInterval = new SimpleInterval(originalInterval.getContig(), Math.max(1, originalInterval.getStart() - fragmentSize), originalInterval.getEnd() + fragmentSize);
        Utils.stream(readsContext.iterator(expandedInterval)).forEach(mate -> {
            final GATKRead read = (mate.isFirstOfPair() ? readTwos : readOnes).get(mate.getName());
            if (read != null) {
                result.put(read, mate);
            }
        });

        return result;
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

    /**
     * Find the 0-based index within a read base array corresponding to a given 1-based position in the reference, along with the cigar operator of
     * the element containing that base.  If the reference coordinate occurs within a deletion, the first index after the deletion is returned.
     * Note that this treats soft-clipped bases as if they align with the reference, which is useful for hard-clipping reads with soft clips.
     *
     * @param alignmentStart        The soft start of the read on the reference
     * @param cigar                 The read's cigar
     * @param refCoord              The target reference coordinate
     * @return                      If the reference coordinate occurs before the read start or after the read end {@code CLIPPING_GOAL_NOT_REACHED};
     *                              if the reference coordinate falls within an alignment block of the read's cigar, the corresponding read coordinate;
     *                              if the reference coordinate falls within a deletion, the first read coordinate after the deletion.  Note: if the last cigar element is
     *                              a deletion (which isn't meaningful), it returns {@code CLIPPING_GOAL_NOT_REACHED}.
     */
    public static Pair<Integer, CigarOperator> getReadIndexForReferenceCoordinate(final int alignmentStart, final Cigar cigar, final int refCoord) {
        if (refCoord < alignmentStart) {
            return new MutablePair<>(READ_INDEX_NOT_FOUND, null);
        }
        int firstReadPosOfElement = 0;              //inclusive
        int firstRefPosOfElement = alignmentStart;  //inclusive
        int lastReadPosOfElement = 0;               //exclusive
        int lastRefPosOfElement = alignmentStart;   //exclusive

        // advance forward through all the cigar elements until we bracket the reference coordinate
        for (final CigarElement element : cigar) {
            final CigarOperator operator = element.getOperator();
            firstReadPosOfElement = lastReadPosOfElement;
            firstRefPosOfElement = lastRefPosOfElement;
            lastReadPosOfElement += operator.consumesReadBases() ? element.getLength() : 0;
            lastRefPosOfElement += operator.consumesReferenceBases() || operator == CigarOperator.S ? element.getLength() : 0;

            if (firstRefPosOfElement <= refCoord && refCoord < lastRefPosOfElement) {   // refCoord falls within this cigar element
                final int readPosAtRefCoord = firstReadPosOfElement + (operator.consumesReadBases() ? ( refCoord - firstRefPosOfElement) : 0);
                return Pair.of(readPosAtRefCoord, operator);
            }
        }
        return new MutablePair<>(READ_INDEX_NOT_FOUND, null);
    }


    /**
     * Returns the index within the read's bases array corresponding to the requested reference coordinate -- or the read coordinate immediately preceding
     * a deletion in which the reference coordinate falls -- along with the cigar operator in which the reference coordinate occurs.
     */
    public static Pair<Integer, CigarOperator> getReadIndexForReferenceCoordinate(final GATKRead read, final int refCoord) {
        return getReadIndexForReferenceCoordinate(read.getSoftStart(), read.getCigar(), refCoord);
    }

    public static Optional<Byte> getReadBaseAtReferenceCoordinate(final GATKRead read, final int refCoord) {
        if (refCoord < read.getStart() || read.getEnd() < refCoord) {
            return Optional.empty();
        }
        final Pair<Integer, CigarOperator> offsetAndOperator = getReadIndexForReferenceCoordinate(read, refCoord);
        return (offsetAndOperator.getLeft() != READ_INDEX_NOT_FOUND && offsetAndOperator.getRight().consumesReadBases()) ?
                Optional.of(read.getBase(offsetAndOperator.getLeft())) : Optional.empty();
    }

    public static Optional<Byte> getReadBaseQualityAtReferenceCoordinate(final GATKRead read, final int refCoord) {
        if (refCoord < read.getStart() || read.getEnd() < refCoord) {
            return Optional.empty();
        }
        final Pair<Integer, CigarOperator> offsetAndOperator = getReadIndexForReferenceCoordinate(read.getSoftStart(), read.getCigar(), refCoord);
        return (offsetAndOperator.getRight() != null && offsetAndOperator.getRight().consumesReadBases()) ?
                Optional.of(read.getBaseQuality(offsetAndOperator.getLeft())) : Optional.empty();
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
     * Creates an "empty", unmapped read with the provided read's read group and mate
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

        emptyRead.setIsUnmapped();
        emptyRead.setMappingQuality(0);
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
        return createCommonSAMWriter(
            (null == outputFile ? null : outputFile.toPath()),
            null == referenceFile ? null : referenceFile.toPath(),
            header,
            preSorted,
            createOutputBamIndex,
            createMD5);
    }

    /**
     * Create a common SAMFileWriter for use with GATK tools.
     *
     * @param outputPath - if this file has a .cram extension then a reference is required. Can not be null.
     * @param referenceFile - the reference source to use. Can not be null if a output file has a .cram extension.
     * @param header - header to be used for the output writer
     * @param preSorted - if true then the records must already be sorted to match the header sort order
     * @param createOutputBamIndex - if true an index will be created for .BAM and .CRAM files
     * @param createMD5 - if true an MD5 file will be created
     *
     * @return SAMFileWriter
     */
    public static SAMFileWriter createCommonSAMWriter(
        final Path outputPath,
        final Path referenceFile,
        final SAMFileHeader header,
        final boolean preSorted,
        boolean createOutputBamIndex,
        final boolean createMD5)
    {
        Utils.nonNull(outputPath);
        Utils.nonNull(header);

        if (createOutputBamIndex && header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            logger.warn("Skipping index file creation for: " +
                outputPath +  ". Index file creation requires reads in coordinate sorted order.");
            createOutputBamIndex = false;
        }

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(createOutputBamIndex).setCreateMd5File(createMD5);
        return ReadUtils.createCommonSAMWriterFromFactory(factory, outputPath, referenceFile, header, preSorted);
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
        return createCommonSAMWriterFromFactory(factory,
            Utils.nonNull(outputFile).toPath(), referenceFile == null ? null : referenceFile.toPath(), header, preSorted);
    }

    /**
     * Create a common SAMFileWriter from a factory for use with GATK tools. Assumes that if the factory has been set
     * to create an index, the header must be set to coordinate sorted.
     *
     * @param outputPath if this file has a .cram extension then a reference is required. Can not be null.
     * @param referenceFile the reference source to use. Can not be null if a output file has a .cram extension.
     * @param header header to be used for the output writer
     * @param preSorted if true then records must already be sorted to match the header sort order
     * @param factory SAMFileWriterFactory factory to use
     * @param openOptions (optional) NIO options specifying how to open the file
     * @return SAMFileWriter
     */
    public static SAMFileWriter createCommonSAMWriterFromFactory(
        final SAMFileWriterFactory factory,
        final Path outputPath,
        final Path referenceFile,
        final SAMFileHeader header,
        final boolean preSorted,
        OpenOption... openOptions)
    {
        Utils.nonNull(outputPath);
        Utils.nonNull(header);

        if (null == referenceFile && outputPath.toString().endsWith(FileExtensions.CRAM)) {
            throw new UserException.MissingReference("A reference file is required for writing CRAM files");
        }

        return factory.makeWriter(header.clone(), preSorted, outputPath, referenceFile);
    }

    /**
     * Validate that a file has CRAM contents by checking that it has a valid CRAM file header
     * (no matter what the extension).
     *
     * @param putativeCRAMPath File to check.
     * @return true if the file has a valid CRAM file header, otherwise false
     */
    public static boolean hasCRAMFileContents(final Path putativeCRAMPath) {
        try (final InputStream fileStream = Files.newInputStream(putativeCRAMPath)) {
            try (final BufferedInputStream bis = new BufferedInputStream(fileStream)) {
                return SamStreams.isCRAMFile(bis);
            }
        }
        catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(e.getMessage());
        }
    }

    /**
     * Validate that a file has CRAM contents by checking that it has a valid CRAM file header
     * (no matter what the extension).
     *
     * @param putativeCRAMFile File to check.
     * @return true if the file has a valid CRAM file header, otherwise false
     */
    public static boolean hasCRAMFileContents(final File putativeCRAMFile) {
        return hasCRAMFileContents(putativeCRAMFile.toPath());
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
        final int adaptorBoundary = read.getAdaptorBoundary();
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
     *
     * @throws IllegalArgumentException if the input {@code read} is {@code null} or does not have any base aligned
     *  against the reference (e.g. is unmapped).
     *
     * @return a number between 0 and the read length-1.
     */
    public static int getFirstAlignedBaseOffset(final GATKRead read) {
        Utils.nonNull(read, "the input read cannot be null");
        if (read.isUnmapped()) {
            throw new IllegalArgumentException("the input read is unmapped and therefore does not have any base aligned");
        } else {
            final List<CigarElement> cigarElements = read.getCigarElements();
            if (cigarElements.isEmpty()) {
                throw new IllegalArgumentException("the input read is mapped yet contains no cigar-elements: " + read.commonToString());
            }
            int result = 0;
            for (final CigarElement ce : cigarElements) {
                final int length = ce.getLength();
                final CigarOperator co = ce.getOperator();
                if (length > 0 && co.isAlignment()) {
                    return result;
                } else if (co.consumesReadBases()) {
                    result += length;
                }
            }
            throw new IllegalArgumentException("the input read cigar does not contain any alignment element");
        }
    }

    /**
     * @param read a GATK read
     * @return true if the read is F2R1, false otherwise
     */
    public static boolean isF2R1(final GATKRead read) {
        return read.isReverseStrand() == read.isFirstOfPair();
    }

    /**
     * @param read a GATK read
     * @return true if the read is F1R2, false otherwise
     */
    public static boolean isF1R2(final GATKRead read) {
        return read.isReverseStrand() != read.isFirstOfPair();
    }

    /**
     * Used to be called isUsableRead()
     **/
    public static boolean readHasReasonableMQ(final GATKRead read){
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

}
