package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Wrapper class for using the PathSeq Bwa aligner class in Spark. Encapsulates closing the index when done.
 */
public final class PSBwaAlignerSpark implements AutoCloseable {

    final PSBwaArgumentCollection bwaArgs;
    private final JavaSparkContext ctx;

    public PSBwaAlignerSpark(final JavaSparkContext ctx, final PSBwaArgumentCollection bwaArgs) {
        this.ctx = ctx;
        this.bwaArgs = bwaArgs;
    }

    public JavaRDD<GATKRead> doBwaAlignment(final JavaRDD<GATKRead> reads,
                                            final boolean pairedAlignment,
                                            final Broadcast<SAMFileHeader> header) {
        final PSBwaArgumentCollection bwaArgsLocal = bwaArgs;
        return reads.mapPartitions(itr -> (new PSBwaAligner(bwaArgsLocal, pairedAlignment)).apply(itr, header.value()));
    }

    //Run this after invoking a Spark action on all RDDs returned from doBwaAlignment()
    public void close() {
        BwaMemIndexCache.closeAllDistributedInstances(ctx);
    }
}
