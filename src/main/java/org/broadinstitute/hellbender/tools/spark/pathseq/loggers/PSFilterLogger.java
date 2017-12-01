package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Interface for filter metrics logging
 */
public interface PSFilterLogger extends AutoCloseable {

    void logPrimaryReads(final JavaRDD<GATKRead> reads);
    void logReadsAfterPrealignedHostFilter(final JavaRDD<GATKRead> reads);
    void logReadsAfterQualityFilter(final JavaRDD<GATKRead> reads);
    void logReadsAfterHostFilter(final JavaRDD<GATKRead> reads);
    void logReadsAfterDeduplication(final JavaRDD<GATKRead> reads);
    void logFinalPairedReads(final JavaRDD<GATKRead> reads);
    @Override
    void close();

}
