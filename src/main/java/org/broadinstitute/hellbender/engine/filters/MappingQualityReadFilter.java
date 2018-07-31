package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep only reads with mapping qualities within a specified range.
 *
 * <p>Note: this filter is not designed to handle the unavailable mapping quality (255).
 * Use MappingQualityAvailableReadFilter to explicitly filter out reads with unavailable quality.</p>
 *
 * @see org.broadinstitute.hellbender.utils.QualityUtils#MAPPING_QUALITY_UNAVAILABLE
 * @see org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.MappingQualityAvailableReadFilter
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads with mapping qualities within a specified range", extraDocs = {ReadFilterLibrary.MappingQualityAvailableReadFilter.class})
public final class MappingQualityReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName=ReadFilterArgumentDefinitions.MINIMUM_MAPPING_QUALITY_NAME, doc = "Minimum mapping quality to keep (inclusive)", optional=true)
    public int minMappingQualityScore = 10;

    @Argument(fullName=ReadFilterArgumentDefinitions.MAXIMUM_MAPPING_QUALITY_NAME, doc = "Maximum mapping quality to keep (inclusive)", optional=true)
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