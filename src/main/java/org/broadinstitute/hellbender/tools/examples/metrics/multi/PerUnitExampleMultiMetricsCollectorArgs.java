package org.broadinstitute.hellbender.tools.examples.metrics.multi;

import java.io.Serializable;

/**
 * Since a given record may be processed by multiple unit collectors (i.e., a sample and a library),
 * the collector args provides an object to hold the result of processing each record so the results
 * of the any calculation can be re-used across all levels.
 *
 * Called by MultiLevelReducibleCollector.makeArg once for each record processed. The result is passed
 * to the sample/library/read-group distributors.
 */
final class PerUnitExampleMultiMetricsCollectorArgs implements Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * The example collector doesn't use this, but could cache some metric/calculation
     * derived from the record.
     */
    public PerUnitExampleMultiMetricsCollectorArgs() {
    }
}
