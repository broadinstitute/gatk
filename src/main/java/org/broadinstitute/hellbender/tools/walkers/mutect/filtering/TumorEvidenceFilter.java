package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.Datum;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.SomaticClusteringModel;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class TumorEvidenceFilter extends Mutect2AlleleFilter {
    @Override
    public ErrorType errorType() { return ErrorType.SEQUENCING; }

    @Override
    protected List<Double> calculateErrorProbabilityForAlleles(final VariantContext vc, final Mutect2FilteringEngine filteringEngine, ReferenceContext referenceContext)
    {
        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final int[] ADs = filteringEngine.sumADsOverSamples(vc, true, false);
        final int totalCount = (int) MathUtils.sum(ADs);
        SomaticClusteringModel model = filteringEngine.getSomaticClusteringModel();

        return IntStream.range(0, tumorLods.length).mapToObj(i ->
                new Datum(tumorLods[i], 0, 0, ADs[i+1], totalCount, SomaticClusteringModel.indelLength(vc, i)))
                .map(model::probabilityOfSequencingError).collect(Collectors.toList());

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
    protected List<String> requiredInfoAnnotations() { return Collections.singletonList(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY); }

}
