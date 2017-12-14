package org.broadinstitute.hellbender.metrics.analysis;

import htsjdk.samtools.metrics.MetricBase;

public final class BaseDistributionByCycleMetrics extends MetricBase {
    public int READ_END;
    public int CYCLE;
    public double PCT_A;
    public double PCT_C;
    public double PCT_G;
    public double PCT_T;
    public double PCT_N;
}
