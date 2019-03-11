package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

public class ReadOrientationFilter extends Mutect2VariantFilter {
    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {

        if (! vc.isSNP()){
            return 0;
        }

        final List<ImmutablePair<Integer, Double>> depthsAndPosteriors = new ArrayList<>();

        vc.getGenotypes().stream().filter(filteringEngine::isTumor)
                .filter(g -> g.hasExtendedAttribute(GATKVCFConstants.ROF_POSTERIOR_KEY) && g.hasExtendedAttribute(GATKVCFConstants.ROF_PRIOR_KEY))
                .forEach(g -> {
                    final double artifactPosterior = GATKProtectedVariantContextUtils.getAttributeAsDouble(g, GATKVCFConstants.ROF_POSTERIOR_KEY, 0.0);
                    final int[] ADs = g.getAD();
                    final int altCount = (int) MathUtils.sum(ADs) - ADs[0];

                    depthsAndPosteriors.add(ImmutablePair.of(altCount, artifactPosterior));
                });

        final double artifactPosterior = weightedMedianPosteriorProbability(depthsAndPosteriors);
        return artifactPosterior;
    }

    @Override
    public String filterName() { return GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME; }

    // the posterior is already annotated in the genotypes.  There's no variant-level posterior.
    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.empty();
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }
}
