package org.broadinstitute.hellbender.tools.variantdb.arrays.tables;

import java.util.Arrays;
import java.util.List;

public class ProbeQcMetricsSchema {

    public static final String PROBE_ID = "probe_id";
    public static final String EXCESS_HET = "excess_het";
    public static final String CALL_RATE = "call_rate";
    public static final String INVARIANT = "invariant";
    
    public static final List<String> PROBE_QC_METRIC_FIELDS = Arrays.asList(PROBE_ID, EXCESS_HET, CALL_RATE, INVARIANT);

}
