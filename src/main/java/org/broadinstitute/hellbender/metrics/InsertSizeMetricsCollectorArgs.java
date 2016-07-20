package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SamPairUtil;

import java.io.Serializable;

// Arguments that need to be calculated once per SAMRecord that are then passed to each PerUnitMetricCollector
// for the given record
final class InsertSizeMetricsCollectorArgs implements Serializable {

    private static final long serialVersionUID = 1L;

    private final int insertSize;
    private final SamPairUtil.PairOrientation po;

    public int getInsertSize() {
        return insertSize;
    }

    public SamPairUtil.PairOrientation getPairOrientation() {
        return po;
    }

    public InsertSizeMetricsCollectorArgs(final int insertSize, final SamPairUtil.PairOrientation po) {
        this.insertSize = insertSize;
        this.po = po;
    }
}
