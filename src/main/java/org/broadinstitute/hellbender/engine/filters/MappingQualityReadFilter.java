package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;

/**
 * Keep reads with high mapping qualities.
 */
public final class MappingQualityReadFilter implements ReadFilter {
    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for calling", optional = true)
    public int minMappingQualtyScore = 10;

    private static final long serialVersionUID = 1L;
    @Override
    public boolean test(final SAMRecord read) {
        return read.getMappingQuality() >= minMappingQualtyScore;
    }
}