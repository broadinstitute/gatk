package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * The state of the allele-fraction model, containing: <p>
 *      1.  the global mean allelic bias <p>
 *      2.  the global variance of the allelic bias <p>
 *      3.  the global outlier probability <p>
 *      4.  minor-allele fractions for each segment <p>
 * <p>
 * See docs/CNVs/CNV-methods.pdf for details.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionState extends ParameterizedState<AlleleFractionParameter> {
    public static final class MinorFractions extends ArrayList<Double> {
        private static final long serialVersionUID = 1029384756L;
        public MinorFractions(final int numSegments) { super(numSegments); }
        public MinorFractions(final List<Double> other) {
            super(new ArrayList<>(other));
        }
    }

    public AlleleFractionState(final double meanBias, final double biasVariance, final double outlierProbability,
                               final MinorFractions minorFractions) {
        super(Arrays.asList(
                new Parameter<>(AlleleFractionParameter.MEAN_BIAS, meanBias),
                new Parameter<>(AlleleFractionParameter.BIAS_VARIANCE, biasVariance),
                new Parameter<>(AlleleFractionParameter.OUTLIER_PROBABILITY, outlierProbability),
                new Parameter<>(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, minorFractions)));
    }

    /**
     * Constructs a single-segment state.
     */
    public AlleleFractionState(final double meanBias, final double biasVariance, final double outlierProbability,
                               final double minorFraction) {
        this(meanBias, biasVariance, outlierProbability, new AlleleFractionState.MinorFractions(Collections.singletonList(minorFraction)));
    }

    public double meanBias() {
        return get(AlleleFractionParameter.MEAN_BIAS, Double.class);
    }

    public double biasVariance() {
        return get(AlleleFractionParameter.BIAS_VARIANCE, Double.class);
    }

    public double outlierProbability() {
        return get(AlleleFractionParameter.OUTLIER_PROBABILITY, Double.class);
    }

    public double segmentMinorFraction(final int segment) {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class).get(segment);
    }

    //copy the state with reference to the SAME minorFractions list data (to save copying) but different value of
    //one of the scalar parameters
    //This is dangerous and minorFractions should not be modified in the copy.
    //
    //The purpose of this is to make an MCMC proposal state for calculating a likelihood with one of the scalar parameters
    //modified (these are unboxed in the getters, so changing these in the copy is safe)
    protected AlleleFractionState shallowCopyWithProposedMeanBias(final double proposedMeanBias) {
        return new AlleleFractionState(proposedMeanBias, biasVariance(), outlierProbability(), minorFractions());
    }

    protected AlleleFractionState shallowCopyWithProposedBiasVariance(final double proposedBiasVariance) {
        return new AlleleFractionState(meanBias(), proposedBiasVariance, outlierProbability(), minorFractions());
    }

    protected AlleleFractionState shallowCopyWithProposedOutlierProbability(final double proposedOutlierProbability) {
        return new AlleleFractionState(meanBias(), biasVariance(), proposedOutlierProbability, minorFractions());
    }

    private MinorFractions minorFractions() {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class);
    }
}
