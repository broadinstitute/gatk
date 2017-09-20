package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Created by valentin on 6/1/17.
 */
public class AlignmentPenalties {

    public static final String UNMAPPED_FRAGMENT_PENALTY_FULL_NAME = "unmappedFragmentPenalty";
    public static final String UNMAPPED_FRAGMENT_PENALTY_SHORT_NAME = "ufp";
    public static final String IMPROPER_PAIR_PENALTY_FULL_NAME = "improperPairPenalty";
    public static final String IMPROPER_PAIR_PENALTY_SHORT_NAME = "ipp";
    public static final String MAXIMUM_TEMPLATE_SCORE_DIFF_FULL_NAME = "maximumTemplateScoreDifference";
    public static final String MAXIMUM_TEMPLATE_SCORE_DIFF_SHORT_NAME = "mtsd";

    public static final double DEFAULT_UNMAPPED_FRAGMENT_PENALTY = 20;
    public static final double DEFAULT_IMPROPER_PAIR_PENALTY = 20;
    public static final double MAXIMUM_LIKELIHOOD_DIFFERENCE_CAP_DEFAULT = 50;

    @Argument(doc = "unmapped fragment penalty",
              fullName = UNMAPPED_FRAGMENT_PENALTY_FULL_NAME,
              shortName = UNMAPPED_FRAGMENT_PENALTY_SHORT_NAME,
              optional = true)
    public double unmappedFragmentPenalty = DEFAULT_UNMAPPED_FRAGMENT_PENALTY;


    @Argument(doc = "improper pair penalty",
              fullName = IMPROPER_PAIR_PENALTY_FULL_NAME,
              shortName = IMPROPER_PAIR_PENALTY_SHORT_NAME,
              optional = true)
    public double improperPairPenalty = DEFAULT_IMPROPER_PAIR_PENALTY;

    @Argument(doc = "maximum fragment score difference ",
             fullName = MAXIMUM_TEMPLATE_SCORE_DIFF_FULL_NAME,
             shortName = MAXIMUM_TEMPLATE_SCORE_DIFF_SHORT_NAME)
    public double maximumTemplateScoreDifference = MAXIMUM_LIKELIHOOD_DIFFERENCE_CAP_DEFAULT;

    public double inversion = -4.0 * Math.log(10);
    public double indelStart = -4.5 * Math.log(10);
    public double indelExtend = Math.log(0.1);
    public double maximumLikelihoodDiffernenceCap;

    public void validate() {
        if (unmappedFragmentPenalty < 0 || Double.isInfinite(unmappedFragmentPenalty) || Double.isNaN(unmappedFragmentPenalty)) {
            throw new CommandLineException.BadArgumentValue(UNMAPPED_FRAGMENT_PENALTY_FULL_NAME, ""+ unmappedFragmentPenalty);
        }
        if (improperPairPenalty < 0 || Double.isInfinite(improperPairPenalty) || Double.isNaN(improperPairPenalty)) {
            throw new CommandLineException.BadArgumentValue(IMPROPER_PAIR_PENALTY_FULL_NAME, "" + improperPairPenalty);
        }
    }

}
