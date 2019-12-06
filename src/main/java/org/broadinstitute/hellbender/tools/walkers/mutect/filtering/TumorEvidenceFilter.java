package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.Datum;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.SomaticClusteringModel;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class TumorEvidenceFilter extends Mutect2AlleleFilter<Integer> {
    @Override
    public ErrorType errorType() { return ErrorType.SEQUENCING; }

    @Override
    protected List<Double> calculateErrorProbabilityForAlleles(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext)
    {
        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final int[] ADs = filteringEngine.sumADsOverSamples(vc, true, false);
        final int totalCount = (int) MathUtils.sum(ADs);
        SomaticClusteringModel model = filteringEngine.getSomaticClusteringModel();

        List<Double> altResults = new ArrayList<>();
        new IndexRange(0, tumorLods.length).forEach(i ->
                altResults.add(model.probabilityOfSequencingError(new Datum(tumorLods[i], 0, 0, ADs[i+1], totalCount, SomaticClusteringModel.indelLength(vc, i)))));

        return altResults;
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
