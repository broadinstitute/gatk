package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.util.ArrayList;
import java.util.List;

/**
 * The collection of those arguments for ApplyBQSR that are not already defined in RecalibrationArgumentCollection.
 */
public class ApplyBQSRUniqueArgumentCollection implements ArgumentCollectionDefinition {
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
     * Static quantized quals are entirely separate from the quantize_qual option which uses dynamic binning.
     * The two types of binning should not be used together.
     */
    @Advanced
    @Argument(fullName="static_quantized_quals", shortName = "SQQ", doc = "Use static quantized quality scores to a given number of levels (with -"+ StandardArgumentDefinitions.BQSR_TABLE_SHORT_NAME+ ")", optional=true, mutex = "quantize_quals")
    public List<Integer> staticQuantizationQuals = new ArrayList<>();

    /**
     * Round down quantized only works with the static_quantized_quals option, and should not be used with
     * the dynamic binning option provided by quantize_quals.  When roundDown = false, rounding is done in
     * probability space to the nearest bin.  When roundDown = true, the value is rounded to the nearest bin
     * that is smaller than the current bin.
     */
    @Advanced
    @Argument(fullName="round_down_quantized", shortName = "RDQ", doc = "Round quals down to nearest quantized qual", optional=true, mutex = "quantize_quals")
    public boolean roundDown = false;

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
        ret.emitOriginalQuals = this.emitOriginalQuals;
        ret.PRESERVE_QSCORES_LESS_THAN = PRESERVE_QSCORES_LESS_THAN;
        ret.globalQScorePrior = this.globalQScorePrior;
        return ret;
    }
}
