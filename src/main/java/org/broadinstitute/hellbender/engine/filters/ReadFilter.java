package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.util.function.Predicate;

public interface ReadFilter extends Predicate<SAMRecord> {

    /*
    Standard read filters
     */
    static ReadFilter UNMAPPED =  read -> !(read.getReadUnmappedFlag() || read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START);
    static ReadFilter NOT_PRIMARY_ALIGNMENT = (samRecord) -> !samRecord.getNotPrimaryAlignmentFlag();
    static ReadFilter DUPLICATE = (samRecord) -> !samRecord.getDuplicateReadFlag();
    static ReadFilter FAILS_VENDOR_QUALITY_CHECK = (samRecord) -> !samRecord.getReadFailsVendorQualityCheckFlag();
    static ReadFilter MAPPING_QUALITY_UNAVAIALBLE = read -> read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    static ReadFilter MAPPING_QUALITY_ZERO = read -> read.getMappingQuality() != 0;

    //HACK: These methods are a hack to get to get the type system to accept compositions of ReadFilters.
    default ReadFilter and(ReadFilter filter ) {
        return Predicate.super.and(filter)::test;
    }

    default ReadFilter or(ReadFilter filter ) {
        return Predicate.super.or(filter)::test;
    }

    default ReadFilter negate(){
        return Predicate.super.negate()::test;
    }
}
