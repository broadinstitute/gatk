package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Keep reads with mapping qualities above a specified threshold.
 */
public final class MappingQualityReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;

    private int minMappingQualityScore = 10;

    public MappingQualityReadFilter( final int minMappingQualityScore ) {
        this.minMappingQualityScore = minMappingQualityScore;
    }

    @Override
    public boolean test( final GATKRead read ) {
        return read.getMappingQuality() >= minMappingQualityScore;
    }
}