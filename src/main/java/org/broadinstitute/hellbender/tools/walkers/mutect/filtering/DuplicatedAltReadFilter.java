package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.UniqueAltReadCount;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

// This filter checks for the case in which PCR-duplicates with unique UMIs (which we assume is caused by false adapter priming)
// amplify the erroneous signal for an alternate allele.
public class DuplicatedAltReadFilter extends HardFilter {
    private final int uniqueAltReadCount;

    public DuplicatedAltReadFilter(final int uniqueAltReadCount) {
        this.uniqueAltReadCount = uniqueAltReadCount;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        return vc.getAttributeAsInt(UniqueAltReadCount.KEY, 1) <= uniqueAltReadCount;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(UniqueAltReadCount.KEY); }
}
