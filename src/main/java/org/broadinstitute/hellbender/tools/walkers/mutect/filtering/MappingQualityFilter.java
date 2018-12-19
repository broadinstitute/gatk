package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.MappingQuality;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class MappingQualityFilter extends HardFilter {
    private final double minMedianMappingQuality;
    private final int longIndelSize;

    public MappingQualityFilter(final double minMedianMappingQuality, final int longIndelSize) {
        this.minMedianMappingQuality = minMedianMappingQuality;
        this.longIndelSize = longIndelSize;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final List<Integer> indelLengths = vc.getIndelLengths();
        final int indelLength = indelLengths == null ? 0 : indelLengths.stream().mapToInt(Math::abs).max().orElseGet(() -> 0);
        final List<Integer> mappingQualityByAllele = vc.getAttributeAsIntList(MappingQuality.KEY, 0);

        // we use the mapping quality annotation of the alt allele in most cases, but for long indels we use the reference
        // annotation.  We have to do this because the indel, even if it maps uniquely, gets a poor mapping quality
        // by virtue of its mismatch.  The reference mapping quality is a decent proxy for the region's mappability.
        return mappingQualityByAllele.get(indelLength < longIndelSize ? 1 : 0) < minMedianMappingQuality;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.MEDIAN_MAPPING_QUALITY_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(MappingQuality.KEY); }
}
