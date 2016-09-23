package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn.AllelicCountTableVerbosity;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Represents an {@link AllelicCount} with posterior probabilities for each site being ref minor, alt minor,
 * or an outlier, according to a model fit by {@link AlleleFractionModeller}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCountWithPhasePosteriors extends AllelicCount {
    private final double refMinorProb;
    private final double altMinorProb;
    private final double outlierProb;

    /**
     * Construct the {@link AllelicCountWithPhasePosteriors} object.  If unnormalized probabilities are passed,
     * they will be normalized.
     * @param count         {@link AllelicCount} of any verbosity (see {@link AllelicCountTableVerbosity})
     * @param refMinorLogProb  log posterior probability of ref-minor phase (can be unnormalized)
     * @param altMinorLogProb  log posterior probability of alt-minor phase (can be unnormalized)
     * @param outlierLogProb   log posterior probability of outlier phase (can be unnormalized)
     */
    public AllelicCountWithPhasePosteriors(final AllelicCount count,
                                           final double refMinorLogProb, final double altMinorLogProb, final double outlierLogProb) {
        super(count);
        ParamUtils.isFinite(refMinorLogProb, "Cannot construct AllelicCountWithPhasePosteriors with non-finite ref-minor probability at: " + count.getInterval());
        ParamUtils.isFinite(altMinorLogProb, "Cannot construct AllelicCountWithPhasePosteriors with non-finite alt-minor probability at: " + count.getInterval());
        ParamUtils.isFinite(outlierLogProb, "Cannot construct AllelicCountWithPhasePosteriors with non-finite outlier probability at: " + count.getInterval());
        final double[] normalizedProbs = MathUtils.normalizeFromLog10ToLinearSpace(
                new double[]{MathUtils.logToLog10(refMinorLogProb), MathUtils.logToLog10(altMinorLogProb), MathUtils.logToLog10(outlierLogProb)});
        refMinorProb = normalizedProbs[0];
        altMinorProb = normalizedProbs[1];
        outlierProb = normalizedProbs[2];
        Utils.validateArg(Double.isFinite(refMinorProb) && Double.isFinite(altMinorProb) && Double.isFinite(outlierProb),
                "Cannot construct AllelicCountWithPhasePosteriors with non-finite probabilities at: " + count.getInterval());
    }

    public double getRefMinorProb() {
        return refMinorProb;
    }

    public double getAltMinorProb() {
        return altMinorProb;
    }

    public double getOutlierProb() {
        return outlierProb;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }

        final AllelicCountWithPhasePosteriors that = (AllelicCountWithPhasePosteriors) o;

        return Double.compare(that.refMinorProb, refMinorProb) == 0 &&
                Double.compare(that.altMinorProb, altMinorProb) == 0 &&
                Double.compare(that.outlierProb, outlierProb) == 0;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        long temp;
        temp = Double.doubleToLongBits(refMinorProb);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(altMinorProb);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(outlierProb);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
