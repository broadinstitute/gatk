package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionSegmenter extends ScalarHMMSegmenter<AllelicCount> {

    private final double outlierProbability = AlleleFractionHMM.DEFAULT_OUTLIER_PROBABILITY;

    /**
     * Initialize the segmenter with its data and panel of normals, giving equal weight to a set of evenly-spaced
     * hidden minor allele fraction values.
     *  @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     * @param acc               The {@link AllelicCountCollection} data attached to this segmenter
     */
    public AlleleFractionSegmenter(final int initialNumStates, final AllelicCountCollection acc) {
        super(acc.getCounts().stream().map(AllelicCount::getInterval).sorted(LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList()),
                acc.getCounts().stream().sorted(LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList()), initialMinorFractions(initialNumStates));
    }

    /**
     * evenly-spaced minor allele fractions going from 0 (inclusive) to 1/2 (exclusive)
     * @param K the initial number of hidden states
     */
    private static List<Double> initialMinorFractions(final int K) {
        ParamUtils.isPositive(K, "must have at least one state");
        return Doubles.asList(GATKProtectedMathUtils.createEvenlySpacedPoints(0.001, 0.5, K));
    }

    public List<ModeledSegment> getModeledSegments() {
        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(positions);
        final List<Pair<SimpleInterval, Double>> segmentation = findSegments();
        return segmentation.stream()
                .map(pair -> new ModeledSegment(pair.getLeft(), tc.targetCount(pair.getLeft()), pair.getRight()))
                .collect(Collectors.toList());
    }

    @Override
    protected ClusteringGenomicHMM<AllelicCount, Double> makeModel() {
        return new AlleleFractionHMM(getStates(), getMemoryLength());
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        // TODO: could relearn outlier probability if desired via pseudocounts
    }

    @Override
    protected double minHiddenStateValue() { return AlleleFractionState.MIN_MINOR_FRACTION; }

    @Override
    protected double maxHiddenStateValue() { return  AlleleFractionState.MAX_MINOR_FRACTION; }

    public double getOutlierProbability() { return outlierProbability; }
}
