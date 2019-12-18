package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;

/**
 * This class stores naming standards in the {@link GermlineCNVCaller}. These constants are also defined in
 * python gcnvkernel package.
 */
public final class GermlineCNVNamingConstants {
    public final static String COPY_NUMBER_POSTERIOR_FILE_NAME = "log_q_c_tc.tsv";
    public final static String COPY_NUMBER_SEGMENTS_FILE_NAME = "copy_number_segments.tsv";
    public final static String BASELINE_COPY_NUMBER_FILE_NAME = "baseline_copy_number_t.tsv";
    public final static String SAMPLE_NAME_TXT_FILE = "sample_name.txt";
    public final static String COPY_NUMBER_TABLE_COLUMN_PREFIX = "COPY_NUMBER_";
    public final static String BASELINE_COPY_NUMBER_TABLE_COLUMN = "BASELINE_COPY_NUMBER";
    public final static String CONTIG_COLUMN = "CONTIG";
    public final static String PLOIDY_COLUMN = "PLOIDY";
    public final static String PLOIDY_GQ_COLUMN = "PLOIDY_GQ";
    public final static String SAMPLE_PREFIX = "SAMPLE_";
    public final static String INTERVAL_LIST_FILE_NAME = "interval_list.tsv";
    public final static String DENOISED_COPY_RATIO_MEAN_FILE_NAME = "mu_denoised_copy_ratio_t.tsv";
    public final static String DEFAULT_GCNV_OUTPUT_COLUMN_PREFIX = "VALUE_";
    public final static String CONTIG_PLOIDY_FILE_NAME = "contig_ploidy.tsv";
}
