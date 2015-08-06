package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;

import java.io.Serializable;

/**
 * The reason for this class' existence is that ApplyWholeBQSRDataflow runs both phases of BQSR.
 * As such, it needs the options for phase 1 and the options for phase 2.
 * The command-line parser errors out if a single class contains two different AgumentCollection classes with the same member.
 * Thus I couldn't use the existing ApplyBQSRArgumentCollection and had to make a variant.
 */
public class ApplyBQSRWithoutMinQScoreArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    /**
     * Turns on the base quantization module. It requires a recalibration report.
     *
     * A value of 0 here means "do not quantize".
     * Any value greater than zero will be used to recalculate the quantization using that many levels.
     * Negative values mean that we should quantize using the recalibration report's quantization level.
     */
    @Argument(fullName="quantize_quals", shortName = "qq", doc = "Quantize quality scores to a given number of levels", optional=true)
    public int quantizationLevels = 0;

    /**
     * Turns off printing of the base insertion and base deletion tags.
     * Only the base substitution qualities will be produced.
     */
    @Argument(fullName="disable_indel_quals", shortName = "DIQ", doc = "Disable printing of base insertion and deletion tags", optional=true)
    public boolean disableIndelQuals = false;

    /**
     * By default, the OQ tag in not emitted. Use this flag to include OQ tags in the output BAM file.
     * Note that this may results in significant file size increase.
     */
    @Argument(fullName="emit_original_quals", shortName = "EOQ", doc = "Emit the OQ tag with the original base qualities", optional=true)
    public boolean emitOriginalQuals = false;

    /**
     * If specified, this value will be used as the prior for all mismatch quality scores instead of the actual reported quality score.
     */
    @Argument(fullName = "globalQScorePrior", shortName = "globalQScorePrior", doc = "Global Qscore Bayesian prior to use for BQSR", optional = true)
    public double globalQScorePrior = -1.0;

    public ApplyBQSRArgumentCollection toApplyBQSRArgumentCollection(int PRESERVE_QSCORES_LESS_THAN) {
        ApplyBQSRArgumentCollection ret = new ApplyBQSRArgumentCollection();
        ret.quantizationLevels = this.quantizationLevels;
        ret.disableIndelQuals = this.disableIndelQuals;
        ret.emitOriginalQuals = this.emitOriginalQuals;
        ret.PRESERVE_QSCORES_LESS_THAN = PRESERVE_QSCORES_LESS_THAN;
        ret.globalQScorePrior = this.globalQScorePrior;
        return ret;
    }
}
