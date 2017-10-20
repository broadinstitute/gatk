package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Dummy filter metrics class that does nothing
 */
public final class PSFilterEmptyLogger implements PSFilterLogger {

    public void logPrimaryReads(final JavaRDD<GATKRead> reads) {}
    public void logReadsAfterPrealignedHostFilter(final JavaRDD<GATKRead> reads) {}
    public void logReadsAfterQualityFilter(final JavaRDD<GATKRead> reads) {}
    public void logReadsAfterHostFilter(final JavaRDD<GATKRead> reads) {}
    public void logReadsAfterDeduplication(final JavaRDD<GATKRead> reads) {}
    public void logFinalPairedReads(final JavaRDD<GATKRead> reads) {}
    public void close() {}

}
