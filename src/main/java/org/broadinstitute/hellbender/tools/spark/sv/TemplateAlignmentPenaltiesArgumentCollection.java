package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.barclay.argparser.Argument;

/**
 * Created by valentin on 6/12/17.
 */
public class TemplateAlignmentPenaltiesArgumentCollection {

    public static final String UNMAPPED_READ_PENALTY_FACTOR_SHORT_NAME = "urpf";

    public static final String UNMAPPED_READ_PENALTY_FACTOR_FULL_NAME = "unmappedReadPenaltyFactor";


    public static final String MAXIMUM_INDEL_PENALTY_SHORT_NAME = "mip";

    public static final String MAXIMUM_INDEL_PENALTY_FULL_NAME = "maximumIndelPenaltyFullName";

    public static final String GAPOPEN_PENALTY_SHORT_NAME = "gop";

    public static final String GAPOPEN_PENALTY_FULL_NAME = "gapOpenPenalty";

    public static final String GAPEXTENSION_PENALTY_SHORT_NAME = "gep";

    public static final String GAPEXTENSION_PENALTY_FULL_NAME = "gapExtensionPenalty";

    public static final String MAXIMUM_MISMATCH_PENALTY_SHORT_NAME = "mmp";

    public static final String MAXIMUM_MISMATCH_PENALTY_FULL_NAME = "maximumMismatchPenalty";

    public static final String IMPROPER_READ_PAIR_PENALTY_FACTOR_SHORT_NAME = "irpf";

    public static final String IMPROPER_READ_PAIR_PENALTY_FACTOR_FULL_NAME = "improperReadPairPenaltyFactor";

    public static final String INVERSION_PENALTY_SHORT_NAME = "ip";

    public static final String INVERSION_PENALTY_FULL_NAME = "inversionPenalty";

    public static final String TRANSLOCATION_PENALTY_SHORT_NAME = "tp";

    public static final String TRANSLOCATION_PENALTY_FULL_NAME = "translocationPenalty";

    public static final double DEFAULT_UNMAPPED_READ_PENALTY_FACTOR = 30.0;
    public static final double DEFAULT_GAPOPEN_PENALTY = 45.0;
    public static final double DEFAULT_GAPEXTENSION_PENALTY = 5.0;

    public static final double DEFAULT_MAXIMUM_INDEL_PENALTY = DEFAULT_GAPOPEN_PENALTY + DEFAULT_GAPEXTENSION_PENALTY * 5;
    public static final double DEFAULT_INVERSION_PENALTY = DEFAULT_MAXIMUM_INDEL_PENALTY + 10.0;
    public static final double DEFAULT_TRANSLOATION_PENALTY = DEFAULT_MAXIMUM_INDEL_PENALTY + 10.0;
    public static final double DEFAULT_MAXIMUM_MISMATCH_PENALTY = 60.0;
    public static final double DEFAULT_IMPROPER_READ_PAIR_PENALTY_FACTOR = 10.0;


    @Argument(doc = "factor added to the worst mapping score amongs mapped reasd to calculate the penalty for reads that don't map at all",
              shortName = UNMAPPED_READ_PENALTY_FACTOR_SHORT_NAME,
              fullName = UNMAPPED_READ_PENALTY_FACTOR_FULL_NAME,
              optional = true)
    public double unmappedReadPenaltyFactor = DEFAULT_UNMAPPED_READ_PENALTY_FACTOR;

    @Argument(doc = "maximum penalty imposable to indel or soft-clips",
              shortName = MAXIMUM_INDEL_PENALTY_SHORT_NAME,
              fullName = MAXIMUM_INDEL_PENALTY_FULL_NAME,
              optional = true)
    public double maximumIndelPenalty = DEFAULT_MAXIMUM_INDEL_PENALTY;

    @Argument(doc = "inversion peanlty imposable to apparent inversion between mates or with reads",
              shortName = INVERSION_PENALTY_SHORT_NAME,
              fullName = INVERSION_PENALTY_FULL_NAME,
              optional = true)
    public double inversionPenalty = DEFAULT_INVERSION_PENALTY;

    @Argument(doc = "translocation penalty imposable to apparent inversion between mates or with reads",
              shortName = TRANSLOCATION_PENALTY_SHORT_NAME,
              fullName = TRANSLOCATION_PENALTY_FULL_NAME,
              optional = true)
    public double translocationPenalty = DEFAULT_TRANSLOATION_PENALTY;

    @Argument(doc = "indel or gap openning penalty",
              shortName = GAPOPEN_PENALTY_SHORT_NAME,
              fullName = GAPOPEN_PENALTY_FULL_NAME,
              optional = true)
    public double gapOpenPenalty = DEFAULT_GAPOPEN_PENALTY;

    @Argument(doc = "indel or gap extension penalty",
              shortName = GAPEXTENSION_PENALTY_SHORT_NAME,
              fullName = GAPEXTENSION_PENALTY_FULL_NAME,
              optional = true)
    public double gapExtensionPenalty = DEFAULT_GAPEXTENSION_PENALTY;

    @Argument(doc = "maximum mismatch penalty impossable",
              shortName = MAXIMUM_MISMATCH_PENALTY_SHORT_NAME,
              fullName = MAXIMUM_MISMATCH_PENALTY_FULL_NAME,
              optional = true)
    public double maximumMismatchPenalty = DEFAULT_MAXIMUM_MISMATCH_PENALTY;

    @Argument(doc = "improper pair penalty factor",
              shortName = IMPROPER_READ_PAIR_PENALTY_FACTOR_SHORT_NAME,
              fullName = IMPROPER_READ_PAIR_PENALTY_FACTOR_FULL_NAME,
              optional = true)
    public double improperPairPenaltyFactor = DEFAULT_IMPROPER_READ_PAIR_PENALTY_FACTOR;
}
