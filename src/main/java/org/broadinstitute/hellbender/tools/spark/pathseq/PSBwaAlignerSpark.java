package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexSingleton;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Wrapper class for using the PathSeq Bwa aligner class in Spark. Encapsulates closing the index when done.
 */
public final class PSBwaAlignerSpark implements AutoCloseable {

    private final JavaSparkContext ctx;
    private final Broadcast<PSBwaArgumentCollection> bwaArgsBroadcast;

    public PSBwaAlignerSpark(final JavaSparkContext ctx, final PSBwaArgumentCollection bwaArgs) {
        this.ctx = ctx;
        this.bwaArgsBroadcast = ctx.broadcast(bwaArgs);
    }

    private static JavaRDD<GATKRead> doBwaAlignmentHelper(final JavaRDD<GATKRead> reads,
                                                          final Broadcast<PSBwaArgumentCollection> bwaArgsBroadcast,
                                                          final boolean pairedAlignment,
                                                          final Broadcast<SAMFileHeader> header) {
        return reads.mapPartitions(itr -> (new PSBwaAligner(bwaArgsBroadcast.value(), pairedAlignment)).apply(itr, header.value()));
    }

    public JavaRDD<GATKRead> doBwaAlignment(final JavaRDD<GATKRead> reads,
                                            final boolean pairedAlignment,
                                            final Broadcast<SAMFileHeader> header) {
        return doBwaAlignmentHelper(reads, bwaArgsBroadcast, pairedAlignment, header);
    }

    //Run this after invoking a Spark action on all RDDs returned from doBwaAlignment()
    public void close() {
        BwaMemIndexSingleton.closeAllDistributedInstances(ctx);
        bwaArgsBroadcast.destroy();
    }
}
