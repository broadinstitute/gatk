package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.HashPartitioner;
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

        //Shuffle reads into partitions by read name hash code, then map each partition to a pair of lists, one
        //  containing the paired reads and the other the unpaired reads
        repartitionedReads = reads.mapToPair(read -> new Tuple2<>(read.getName(), read))
                .partitionBy(new HashPartitioner(reads.getNumPartitions()))
                .map(Tuple2::_2)
                .mapPartitions(iter -> mapPartitionsToPairedAndUnpairedLists(iter, readsPerPartitionGuess));
        repartitionedReads.cache();
        isCached = true;
    }

    /**
     * Maps each partition to a Tuple of two Lists, the first containing the paired reads, the second containing unpaired
     */
    private static Iterator<Tuple2<List<GATKRead>, List<GATKRead>>> mapPartitionsToPairedAndUnpairedLists(final Iterator<GATKRead> iter, final int readsPerPartitionGuess) {
        //Find the paired and unpaired reads by scanning the partition for repeated names
        final List<GATKRead> pairedReadsList = new ArrayList<>(readsPerPartitionGuess);
        final Map<String, GATKRead> unpairedReads = new HashMap<>(readsPerPartitionGuess);
        while (iter.hasNext()) {
            final GATKRead read = iter.next();
            final String readName = read.getName();
            //If read's mate is already in unpairedReads then we have a pair, which gets added to the ordered List
            if (unpairedReads.containsKey(readName)) {
                pairedReadsList.add(read);
                pairedReadsList.add(unpairedReads.remove(readName));
            } else {
                unpairedReads.put(readName, read);
            }
        }
        //Get the unpaired reads out of the hashmap
        final List<GATKRead> unpairedReadsList = new ArrayList<>(unpairedReads.values());

        //Minimize unpairedReads memory footprint (don't rely on readsPerPartitionGuess)
        final List<GATKRead> pairedReadsListResized = new ArrayList<>(pairedReadsList.size());
        pairedReadsListResized.addAll(pairedReadsList);

        //Wrap and return the result
        return Collections.singletonList(new Tuple2<>(pairedReadsListResized, unpairedReadsList)).iterator();
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
