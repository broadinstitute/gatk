package org.broadinstitute.hellbender.tools.walkers.bqsr;

/**
 * Collection of datafile names for BQSR tests.
 */
public final class BQSRTestData {
    public static final String EXPECTED_WGS_B37_CH20_1READ_RECAL = "expected.overlappingRead.recal.txt";
    public static final String EXPECTED_WGS_B37_CH20_1READ_NOREFBASES_RECAL = "expected.overlappingRead.noRefBases.recal.txt.gz";
    public static final String EXPECTED_WGS_B37_CH20_1M_1M1K_RECAL = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recal.txt";

    //Note: this recal table was created using GATK4 because this functionality does not exist in GATK3
    public static final String EXPECTED_WGS_B37_CH20_1M_1M1K_NOINDEL_NOBAQ_RECAL = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.noIndels.noBAQ.recal.txt";

    public static final String EXPECTED_WGS_B37_CH20_1M_1M1K_INDELS_CONTEXT_SIZE_4_RECAL = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.indels_context_size_4.recal.txt";
    public static final String EXPECTED_WGS_B37_CH20_1M_1M1K_LOW_QUALITY_TAIL_5_RECAL = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.low_quality_tail_5.recal.txt";
    public static final String EXPECTED_WGS_B37_CH20_1M_1M1K_QUANTIZING_LEVELS_6_RECAL = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.quantizing_levels_6.recal.txt";
    public static final String EXPECTED_WGS_B37_CH20_1M_1M1K_MISMATCHES_CONTEXT_SIZE_4_RECAL = "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.mismatches_context_size_4.recal.txt";
}
