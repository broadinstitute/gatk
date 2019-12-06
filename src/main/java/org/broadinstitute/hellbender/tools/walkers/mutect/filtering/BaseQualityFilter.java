package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class BaseQualityFilter extends HardAlleleFilter<Integer> {
    private final double minMedianBaseQuality;

    public BaseQualityFilter(final double minMedianBaseQuality) {
        this.minMedianBaseQuality = minMedianBaseQuality;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        List<Integer> baseQualityByAllele = vc.getAttributeAsIntList(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY, 0);
        baseQualityByAllele.remove(0); // get rid of ref
        return baseQualityByAllele.stream().map(qual -> qual < minMedianBaseQuality).collect(Collectors.toList());
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY); }
}
