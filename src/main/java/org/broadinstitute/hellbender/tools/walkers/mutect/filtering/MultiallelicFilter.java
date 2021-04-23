package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class MultiallelicFilter extends HardFilter {
    private final int numAltAllelesThreshold;
    private static final double MULTIALLELIC_LOD_THRESHOLD = 5.0;

    public MultiallelicFilter(final int numAltAllelesThreshold) {
        this.numAltAllelesThreshold = numAltAllelesThreshold;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final long numPassingAltAlleles = Arrays.stream(tumorLods)
                .filter(lod -> lod > MULTIALLELIC_LOD_THRESHOLD)
                .count();

        return numPassingAltAlleles > numAltAllelesThreshold;
    }

    public String filterName() {
        return GATKVCFConstants.MULTIALLELIC_FILTER_NAME;
    }

    protected List<String> requiredInfoAnnotations() { return Collections.singletonList(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY); }
}
