package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep reads with mapping qualities within a specified range.
 *
 * Note: this filter does not handle specially the unavailable mapping quality ({@link org.broadinstitute.hellbender.utils.QualityUtils#MAPPING_QUALITY_UNAVAILABLE}).
 * Use {@link org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.MappingQualityAvailableReadFilter} to explicitly filter out reads with unavailable quality.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY)
public final class MappingQualityReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName="minimumMappingQuality", shortName="minimumMappingQuality", doc = "Minimum mapping quality to keep (inclusive)", optional=true)
    public int minMappingQualityScore = 10;

    @Argument(fullName="maximumMappingQuality", shortName="maximumMappingQuality", doc = "Maximum mapping quality to keep (inclusive)", optional=true)
    public Integer maxMappingQualityScore = null;

    // Command line parser requires a no-arg constructor
    public MappingQualityReadFilter() {}

    public MappingQualityReadFilter( final int minMappingQualityScore ) {
        this.minMappingQualityScore = minMappingQualityScore;
    }

    public MappingQualityReadFilter( final int minMappingQualityScore, final Integer maxMappingQualityScore) {
        this.minMappingQualityScore = minMappingQualityScore;
        this.maxMappingQualityScore = maxMappingQualityScore;
    }

    @Override
    public boolean test( final GATKRead read ) {
        final int mq = read.getMappingQuality();
        return  mq >= minMappingQualityScore
                && (maxMappingQualityScore == null || mq <= maxMappingQualityScore);
    }
}