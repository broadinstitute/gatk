package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class ReadPositionFilter extends HardAlleleFilter<Integer> {
    private final double minMedianReadPosition;

    public ReadPositionFilter(final double minMedianReadPosition) {
        this.minMedianReadPosition = minMedianReadPosition;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        // MPOS doesn't have data for ref allele
        final List<Integer> readPositionByAllele = vc.getAttributeAsIntList(GATKVCFConstants.MEDIAN_READ_POSITON_KEY, 0);
        return readPositionByAllele.subList(0, readPositionByAllele.size()).stream()
                // a negative value is possible due to a bug: https://github.com/broadinstitute/gatk/issues/5492
                .map(readPos -> readPos > -1 && readPos < minMedianReadPosition).collect(Collectors.toList());
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.READ_POSITION_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.MEDIAN_READ_POSITON_KEY); }
}
