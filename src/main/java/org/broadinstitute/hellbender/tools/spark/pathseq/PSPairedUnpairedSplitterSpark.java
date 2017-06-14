package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;

/**
 * Class for separating paired and unpaired reads in an RDD
 */
public final class PSPairedUnpairedSplitterSpark {

    protected final Logger logger = LogManager.getLogger(this.getClass());
    private final JavaRDD<Tuple2<List<GATKRead>, List<GATKRead>>> repartitionedReads;
    private boolean isCached;

    /**
     * Gets RDDs of the paired and unpaired reads
     */
    public PSPairedUnpairedSplitterSpark(final JavaRDD<GATKRead> reads, final int readsPerPartitionGuess) {

        //Repartition reads then map each partition to a pair of lists, one containing the paired reads and the
        // other the unpaired reads
        repartitionedReads = PSFilter.repartitionPairedReads(reads)
                .mapPartitions(iter -> mapPartitionsToPairedAndUnpairedLists(iter, readsPerPartitionGuess));
        repartitionedReads.cache();
        isCached = true;
    }

    /**
     * Wrapper for getPairedAndUnpairedLists()
     */
    private static Iterator<Tuple2<List<GATKRead>, List<GATKRead>>> mapPartitionsToPairedAndUnpairedLists(final Iterator<GATKRead> iter, final int readsPerPartitionGuess) {
        //Wrap and return the result
        return Collections.singletonList(PSFilter.getPairedAndUnpairedLists(iter,readsPerPartitionGuess)).iterator();
    }

    public JavaRDD<GATKRead> getPairedReads() {
        if (!isCached) {
            logger.warn("Getting paired reads after call to close(). Performance may be reduced.");
        }
        return repartitionedReads.flatMap(tuple -> tuple._1.iterator());
    }

    public JavaRDD<GATKRead> getUnpairedReads() {
        if (!isCached) {
            logger.warn("Getting unpaired reads after call to close(). Performance may be reduced.");
        }
        return repartitionedReads.flatMap(tuple -> tuple._2.iterator());
    }

    /**
     * Call only after Spark is done with all closures involving the paired/unpaired RDDs,
     * i.e. an action like collect() has been called.
     */
    public void close() {
        repartitionedReads.unpersist();
        isCached = false;
    }
}
