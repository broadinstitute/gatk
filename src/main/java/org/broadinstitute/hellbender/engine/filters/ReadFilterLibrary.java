package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;

/**
 * standard ReadFilters
 */
public final class ReadFilterLibrary {
    public static final ReadFilter ALLOW_ALL_READS = read -> true;

    public static final ReadFilter MAPPED =  read -> !(read.getReadUnmappedFlag() || read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START);
    public static final ReadFilter PRIMARY_ALIGNMENT = read -> !read.getNotPrimaryAlignmentFlag();
    public static final ReadFilter NOT_DUPLICATE = read -> !read.getDuplicateReadFlag();
    public static final ReadFilter PASSES_VENDOR_QUALITY_CHECK = read -> !read.getReadFailsVendorQualityCheckFlag();
    public static final ReadFilter MAPPING_QUALITY_AVAILABLE = read -> read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    public static final ReadFilter MAPPING_QUALITY_NOT_ZERO = read -> read.getMappingQuality() != 0;

    public static final ReadFilter VALID_ALIGNMENT_START = ReadUtils::hasValidAlignmentStart;
    public static final ReadFilter VALID_ALIGNMENT_END = ReadUtils::hasValidAlignmentEnd;
    public static final ReadFilter ALIGNMENT_AGREES_WITH_HEADER = read -> ReadUtils.alignmentAgreesWithHeader(read.getHeader(), read);
    public static final ReadFilter HAS_READ_GROUP = ReadUtils::hasReadGroup;
    public static final ReadFilter HAS_MATCHING_BASES_AND_QUALS = ReadUtils::hasMatchingBasesAndQuals;
    public static final ReadFilter CIGAR_AGREES_WITH_ALIGNEMENT = ReadUtils::cigarAgreesWithAlignement;
    public static final ReadFilter SEQ_IS_STORED = ReadUtils::seqIsStored;
    public static final ReadFilter CIGAR_IS_SUPPORTED = ReadUtils::cigarIsSupported;

    public static final ReadFilter WELLFORMED =
            VALID_ALIGNMENT_START
            .and(VALID_ALIGNMENT_END)
            .and(ALIGNMENT_AGREES_WITH_HEADER)
            .and(HAS_READ_GROUP)
            .and(HAS_MATCHING_BASES_AND_QUALS)
            .and(CIGAR_AGREES_WITH_ALIGNEMENT)
            .and(SEQ_IS_STORED)
            .and(CIGAR_IS_SUPPORTED);

}
