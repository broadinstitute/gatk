package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;

/**
 * This class stores naming standards in the {@link GermlineCNVCaller}. These constants are also defined in
 * python gcnvkernel package.
 */
public final class GermlineCNVNamingConstants {

    public final static String COPY_NUMBER_POSTERIOR_FILE_NAME = "log_q_c_tc.tsv";

    public final static String COPY_NUMBER_TABLE_COLUMN_PREFIX = "COPY_NUMBER_";

    public final static String INTERVAL_LIST_FILE_NAME = "interval_list.tsv";
}
