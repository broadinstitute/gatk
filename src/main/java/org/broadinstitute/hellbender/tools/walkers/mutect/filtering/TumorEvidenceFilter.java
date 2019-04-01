package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.Datum;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.SomaticClusteringModel;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

public class TumorEvidenceFilter extends Mutect2VariantFilter {
    @Override
    public ErrorType errorType() { return ErrorType.SEQUENCING; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext) {
        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final int[] ADs = filteringEngine.sumADsOverSamples(vc, true, false);
        final int maxIndex = MathUtils.maxElementIndex(tumorLods);
        final int altCount = ADs[maxIndex + 1];
        final int totalCount = (int) MathUtils.sum(ADs);

        return filteringEngine.getSomaticClusteringModel()
                .probabilityOfSequencingError(new Datum(tumorLods[maxIndex], 0, 0, altCount, totalCount, SomaticClusteringModel.indelLength(vc, maxIndex)));
    }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.SEQUENCING_QUAL_KEY);
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.TUMOR_EVIDENCE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY); }

}
