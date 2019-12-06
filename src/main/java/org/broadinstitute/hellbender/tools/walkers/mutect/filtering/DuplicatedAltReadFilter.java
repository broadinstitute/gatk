package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.UniqueAltReadCount;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

// This filter checks for the case in which PCR-duplicates with unique UMIs (which we assume is caused by false adapter priming)
// amplify the erroneous signal for an alternate allele.
public class DuplicatedAltReadFilter extends HardAlleleFilter<Integer> {
    private final int uniqueAltReadCount;

    public DuplicatedAltReadFilter(final int uniqueAltReadCount) {
        this.uniqueAltReadCount = uniqueAltReadCount;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public List<Boolean> areAllelesArtifacts(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        return Collections.singletonList(vc.getAttributeAsInt(UniqueAltReadCount.KEY, 1) <= uniqueAltReadCount);
    }

   @Override
    public String filterName() {
        return GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(UniqueAltReadCount.KEY); }
}
