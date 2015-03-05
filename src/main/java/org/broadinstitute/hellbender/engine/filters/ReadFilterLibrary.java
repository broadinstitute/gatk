package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.QualityUtils;

/**
 * standard ReadFilters
 */
public final class ReadFilterLibrary {
    public static ReadFilter MAPPED =  read -> !(read.getReadUnmappedFlag() || read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START);
    public static ReadFilter PRIMARY_ALIGNMENT = read -> !read.getNotPrimaryAlignmentFlag();
    public static ReadFilter NOT_DUPLICATE = read -> !read.getDuplicateReadFlag();
    public static ReadFilter PASSES_VENDOR_QUALITY_CHECK = read -> !read.getReadFailsVendorQualityCheckFlag();
    public static ReadFilter MAPPING_QUALITY_AVAILABLE = read -> read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    public static ReadFilter MAPPING_QUALITY_NOT_ZERO = read -> read.getMappingQuality() != 0;
    public static ReadFilter WELLFORMED = WellformedReads::isWellformed;
    public static ReadFilter ALLOW_ALL_READS = read -> true;
}
