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
    public static final double MIN_MINOR_FRACTION = 0.0;   //by definition!
    public static final double MAX_MINOR_FRACTION = 0.5;   //by definition!

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

    public AlleleFractionGlobalParameters globalParameters() {
        return new AlleleFractionGlobalParameters(meanBias(), biasVariance(), outlierProbability());
    }

    public MinorFractions minorFractions() {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class);
    }
}
