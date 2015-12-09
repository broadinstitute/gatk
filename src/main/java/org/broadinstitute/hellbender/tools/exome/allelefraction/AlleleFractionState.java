package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.utils.mcmc.AbstractParameterizedState;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The state of the allele fraction model, containing: <p>
 *      1.  minor allele fractions for each segment <p>
 *      2.  a global outlier probability <p>
 *      3.  the mean allelic bias <p>
 *      4.  the rate (mean / variance) of the allelic bias <p>
 * <p>
 * See docs/AllelicCapSeg/ACS-methods.pdf for details.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionState extends AbstractParameterizedState {
    @SuppressWarnings("serial")
    public static final class MinorFractions extends ArrayList<Double> {
        public MinorFractions() { super(); }
        public MinorFractions(final List<Double> other) {
            super(new ArrayList<>(other));
        }
    }

    public static final String MEAN_BIAS_NAME = "mu";
    public static final String BIAS_VARIANCE_NAME = "var";
    public static final String P_OUTLIER_NAME = "pi";
    public static final String MINOR_FRACTIONS_NAME = "f";

    @Override
    protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
        return stateClass.cast(new AlleleFractionState(meanBias(), biasVariance(), outlierProbability(),
                new MinorFractions(minorFractions())));
    }

    //copy the state with reference to the SAME minorFractions list data (to save copying) but different value of
    //one of the scalar parameters
    //This is dangerous and minorFractions should not be modified in the copy.
    //
    //The purpose of this is to make an MCMC proposal state for calculating a likelihood with one of the scalar parameters
    //modified (these are unboxed in the getters, so changing these in the copy is safe)
    public AlleleFractionState shallowCopyWithProposedMeanBias(final double proposedMeanBias) {
        return new AlleleFractionState(proposedMeanBias, biasVariance(), outlierProbability(), minorFractions());
    }

    public AlleleFractionState shallowCopyWithProposedBiasVariance(final double proposedBiasVariance) {
        return new AlleleFractionState(meanBias(), proposedBiasVariance, outlierProbability(), minorFractions());
    }

    public AlleleFractionState shallowCopyWithProposedOutlierProbability(final double proposedOutlierProbability) {
        return new AlleleFractionState(meanBias(), biasVariance(), proposedOutlierProbability, minorFractions());
    }

    public AlleleFractionState(final double biasMean, final double biasVariance, final double outlierProbability,
            final MinorFractions minorFractions) {
        super(Arrays.asList(
                new Parameter<>(MEAN_BIAS_NAME, biasMean),
                new Parameter<>(BIAS_VARIANCE_NAME, biasVariance),
                new Parameter<>(P_OUTLIER_NAME, outlierProbability),
                new Parameter<>(MINOR_FRACTIONS_NAME, minorFractions)));
    }

    public double meanBias() {
        return get(MEAN_BIAS_NAME, Double.class);
    }

    public double biasVariance() {
        return get(BIAS_VARIANCE_NAME, Double.class);
    }

    public double outlierProbability() {
        return get(P_OUTLIER_NAME, Double.class);
    }

    public MinorFractions minorFractions() {
        return get(MINOR_FRACTIONS_NAME, MinorFractions.class);
    }

    public double minorFractionInSegment(final int segment) {
        return get(MINOR_FRACTIONS_NAME, MinorFractions.class).get(segment);
    }
}
