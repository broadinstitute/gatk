package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Interface for score metrics logging
 */
public interface PSScoreLogger {

    void logReadCounts(final JavaRDD<GATKRead> reads);
    void close();

}
