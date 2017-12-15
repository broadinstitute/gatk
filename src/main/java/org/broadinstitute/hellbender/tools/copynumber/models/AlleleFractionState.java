package org.broadinstitute.hellbender.tools.copynumber.models;

import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The state of the allele-fraction model, containing: <p>
 *      1.  the global mean reference bias <p>
 *      2.  the global variance of the reference bias <p>
 *      3.  the global outlier probability <p>
 *      4.  minor-allele fractions for each segment <p>
 * <p>
 * See docs/CNVs/CNV-methods.pdf for details.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionState extends ParameterizedState<AlleleFractionParameter> {
    static final class MinorFractions extends ArrayList<Double> {
        private static final long serialVersionUID = 1029384756L;

        MinorFractions(final int numSegments) {
            super(numSegments);
        }

        MinorFractions(final List<Double> minorFractions) {
            super(new ArrayList<>(minorFractions));
        }
    }

    AlleleFractionState(final double meanBias,
                        final double biasVariance,
                        final double outlierProbability,
                        final MinorFractions minorFractions) {
        super(Arrays.asList(
                new Parameter<>(AlleleFractionParameter.MEAN_BIAS, meanBias),
                new Parameter<>(AlleleFractionParameter.BIAS_VARIANCE, biasVariance),
                new Parameter<>(AlleleFractionParameter.OUTLIER_PROBABILITY, outlierProbability),
                new Parameter<>(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, minorFractions)));
    }

    double meanBias() {
        return get(AlleleFractionParameter.MEAN_BIAS, Double.class);
    }

    double biasVariance() {
        return get(AlleleFractionParameter.BIAS_VARIANCE, Double.class);
    }

    double outlierProbability() {
        return get(AlleleFractionParameter.OUTLIER_PROBABILITY, Double.class);
    }

    double segmentMinorFraction(final int segment) {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class).get(segment);
    }

    AlleleFractionGlobalParameters globalParameters() {
        return new AlleleFractionGlobalParameters(meanBias(), biasVariance(), outlierProbability());
    }

    MinorFractions minorFractions() {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class);
    }
}