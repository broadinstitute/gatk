package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Dummy filter metrics class that does nothing
 */
public final class PSFilterEmptyLogger implements PSFilterLogger {

    @Override
    public void logPrimaryReads(final JavaRDD<GATKRead> reads) {}

    @Override
    public void logReadsAfterPrealignedHostFilter(final JavaRDD<GATKRead> reads) {}

    @Override
    public void logReadsAfterQualityFilter(final JavaRDD<GATKRead> reads) {}

    @Override
    public void logReadsAfterHostFilter(final JavaRDD<GATKRead> reads) {}

    @Override
    public void logReadsAfterDeduplication(final JavaRDD<GATKRead> reads) {}

    @Override
    public void logFinalPairedReads(final JavaRDD<GATKRead> reads) {}

    @Override
    public void close() {}

}
