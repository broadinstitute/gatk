package org.broadinstitute.hellbender.tools.copynumber.arguments;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class SomaticGenotypingArgumentCollection implements Serializable {
    public static final long serialVersionUID = 1L;

    //het genotyping argument names
    public static final String MINIMUM_TOTAL_ALLELE_COUNT_CASE_LONG_NAME = "minimum-total-allele-count-case";
    public static final String MINIMUM_TOTAL_ALLELE_COUNT_NORMAL_LONG_NAME = "minimum-total-allele-count-normal";
    public static final String GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD_LONG_NAME = "genotyping-homozygous-log-ratio-threshold";
    public static final String GENOTYPING_BASE_ERROR_RATE_LONG_NAME = "genotyping-base-error-rate";

    @Argument(
            doc = "Minimum total count for filtering allelic counts in the case sample, if available.  " +
                    "The default value of zero is appropriate for matched-normal mode; " +
                    "increase to an appropriate value for case-only mode.",
            fullName = MINIMUM_TOTAL_ALLELE_COUNT_CASE_LONG_NAME,
            minValue = 0,
            optional = true
    )
    public int minTotalAlleleCountCase = 0;

    @Argument(
            doc = "Minimum total count for filtering allelic counts in the matched-normal sample, if available.",
            fullName = MINIMUM_TOTAL_ALLELE_COUNT_NORMAL_LONG_NAME,
            minValue = 0,
            optional = true
    )
    public int minTotalAlleleCountNormal = 30;

    @Argument(
            doc = "Log-ratio threshold for genotyping and filtering homozygous allelic counts, if available.  " +
                    "Increasing this value will increase the number of sites assumed to be heterozygous for modeling.",
            fullName = GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD_LONG_NAME,
            optional = true
    )
    public double genotypingHomozygousLogRatioThreshold = -10.;

    @Argument(
            doc = "Maximum base-error rate for genotyping and filtering homozygous allelic counts, if available.  " +
                    "The likelihood for an allelic count to be generated from a homozygous site will be integrated " +
                    "from zero base-error rate up to this value.  Decreasing this value will increase " +
                    "the number of sites assumed to be heterozygous for modeling.",
            fullName = GENOTYPING_BASE_ERROR_RATE_LONG_NAME,
            minValue = 0.,
            maxValue = 1.,
            optional = true
    )
    public double genotypingBaseErrorRate = 5E-2;
}
