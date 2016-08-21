package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Stores prior information for the {@link TumorHeterogeneity} model of a mixture of subclones with copy-number variation..
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityPriorCollection {
    private final PloidyState normalPloidyState;
    private final PloidyStatePrior ploidyStatePrior;

    private final HyperparameterValues concentrationPriorHyperparameterValues;
    private final HyperparameterValues copyRatioNormalizationPriorHyperparameterValues;
    private final HyperparameterValues copyRatioNoiseConstantPriorHyperparameterValues;

    private final double ploidyMismatchPenalty;
    private final double subcloneVariancePenalty;

    public TumorHeterogeneityPriorCollection(final PloidyState normalPloidyState,
                                             final PloidyStatePrior ploidyStatePrior,
                                             final double concentrationPriorAlpha,
                                             final double concentrationPriorBeta,
                                             final double copyRatioNormalizationPriorAlpha,
                                             final double copyRatioNormalizationPriorBeta,
                                             final double copyRatioNoiseConstantPriorAlpha,
                                             final double copyRatioNoiseConstantPriorBeta,
                                             final double ploidyMismatchPenalty,
                                             final double subcloneVariancePenalty) {
        Utils.nonNull(normalPloidyState);
        Utils.nonNull(ploidyStatePrior);
        Utils.validateArg(ploidyStatePrior.ploidyStates().contains(normalPloidyState),
                "Ploidy-state prior must contain normal ploidy state.");
        Utils.validateArg(ploidyMismatchPenalty >= 0, "Ploidy-mismatch penalty must be non-negative.");
        Utils.validateArg(subcloneVariancePenalty >= 0, "Subclone-variance penalty must be non-negative.");
        this.normalPloidyState = normalPloidyState;
        this.ploidyStatePrior = ploidyStatePrior;
        concentrationPriorHyperparameterValues = new HyperparameterValues(concentrationPriorAlpha, concentrationPriorBeta);
        copyRatioNormalizationPriorHyperparameterValues = new HyperparameterValues(copyRatioNormalizationPriorAlpha, copyRatioNormalizationPriorBeta);
        copyRatioNoiseConstantPriorHyperparameterValues = new HyperparameterValues(copyRatioNoiseConstantPriorAlpha, copyRatioNoiseConstantPriorBeta);
        this.ploidyMismatchPenalty = ploidyMismatchPenalty;
        this.subcloneVariancePenalty = subcloneVariancePenalty;
    }

    public PloidyState normalPloidyState() {
        return normalPloidyState;
    }

    public PloidyStatePrior ploidyStatePrior() {
        return ploidyStatePrior;
    }

    public HyperparameterValues concentrationPriorHyperparameterValues() {
        return concentrationPriorHyperparameterValues;
    }

    public HyperparameterValues copyRatioNormalizationPriorHyperparameterValues() {
        return copyRatioNormalizationPriorHyperparameterValues;
    }

    public HyperparameterValues copyRatioNoiseConstantPriorHyperparameterValues() {
        return copyRatioNoiseConstantPriorHyperparameterValues;
    }

    public double ploidyMismatchPenalty() {
        return ploidyMismatchPenalty;
    }

    public double subcloneVariancePenalty() {
        return subcloneVariancePenalty;
    }

    public static class HyperparameterValues {
        private final double alpha;
        private final double beta;

        public HyperparameterValues(final double alpha, final double beta) {
            Utils.validateArg(alpha > 0, "Hyperparameter must be positive.");
            Utils.validateArg(beta > 0, "Hyperparameter must be positive.");
            this.alpha = alpha;
            this.beta = beta;
        }

        public double getAlpha() {
            return alpha;
        }

        public double getBeta() {
            return beta;
        }
    }
}
