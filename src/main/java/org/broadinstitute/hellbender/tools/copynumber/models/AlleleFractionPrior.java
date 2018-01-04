package org.broadinstitute.hellbender.tools.copynumber.models;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents priors for the allele-fraction model.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionPrior {
    private final double minorAlleleFractionPriorAlpha;

    public AlleleFractionPrior(final double minorAlleleFractionPriorAlpha) {
        Utils.validateArg(minorAlleleFractionPriorAlpha >= 1,
                "Alpha hyperparameter for the 4-parameter beta-distribution prior on " +
                        "segment minor-allele fraction must be greater than or equal to one.");
        this.minorAlleleFractionPriorAlpha = minorAlleleFractionPriorAlpha;
    }

    double getMinorAlleleFractionPriorAlpha() {
        return minorAlleleFractionPriorAlpha;
    }
}
