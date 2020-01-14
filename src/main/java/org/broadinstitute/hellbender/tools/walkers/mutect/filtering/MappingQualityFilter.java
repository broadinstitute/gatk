package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

// TODO this class doens't use the generic method - how to simplify?
public class MappingQualityFilter extends HardAlleleFilter<Integer> {
    private final double minMedianMappingQuality;
    private final int longIndelSize;

    public MappingQualityFilter(final double minMedianMappingQuality, final int longIndelSize) {
        this.minMedianMappingQuality = minMedianMappingQuality;
        this.longIndelSize = longIndelSize;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        final List<Integer> indelLengths = vc.getIndelLengths(); // alts only
        List<Integer> mappingQualityByAllele = vc.getAttributeAsIntList(GATKVCFConstants.MEDIAN_MAPPING_QUALITY_KEY, 0);

        // we use the mapping quality annotation of the alt allele in most cases, but for long indels we use the reference
        // annotation.  We have to do this because the indel, even if it maps uniquely, gets a poor mapping quality
        // by virtue of its mismatch.  The reference mapping quality is a decent proxy for the region's mappability.
        int refQual = mappingQualityByAllele.remove(0); // get the ref value and convert list to alts only
        new IndexRange(0, mappingQualityByAllele.size()).forEach(i -> {
            if (indelLengths != null && indelLengths.get(i) >= longIndelSize) {
                mappingQualityByAllele.set(i, refQual);
            }
        });
        return mappingQualityByAllele.stream().map(qual -> qual < minMedianMappingQuality).collect(Collectors.toList());
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.MEDIAN_MAPPING_QUALITY_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.MEDIAN_MAPPING_QUALITY_KEY); }
}
