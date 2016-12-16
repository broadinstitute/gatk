package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.spark.sv.SVUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

public class PSFilterUtils {

    /**
     * Preferentially filters unpaired reads, when possible. Assumes reads have pairedness flags set properly
     */
    protected static JavaRDD<GATKRead> filterDuplicateSequences(final JavaRDD<GATKRead> reads) {
        return reads.mapToPair(PSFilterUtils::canonicalizeRead)
                .groupByKey()
                .values()
                .map(iter -> {
                    for (final GATKRead read : iter) {
                        if (!read.isPaired()) return read;
                    }
                    return iter.iterator().next();
                });
    }

    /**
     * Pairs reads with canonical Long ID using the lesser 64-bit hash of the sequence and its reverse complement
     */
    protected static Tuple2<Long, GATKRead> canonicalizeRead(final GATKRead read) {
        final byte[] bases = read.getBases();
        final long hashForward = fnvByteArray64(1099511628211L, bases);
        SequenceUtil.reverseComplement(bases);
        final long hashReverse = fnvByteArray64(1099511628211L, bases);
        return new Tuple2<>(Math.min(hashForward, hashReverse), read);
    }

    /**
     * 64-bit FNV-1a hash for byte arrays
     */
    private static long fnvByteArray64(long start, final byte[] toHash) {
        for (int i = 0; i < toHash.length; i += 8) {
            long val = 0;
            for (int j = 0; i + j < toHash.length; j++) {
                val = (val | toHash[j]) << 8;
            }
            start = SVUtils.fnvLong64(start, val);
        }
        return start;
    }

    /**
     * Sets proper pairedness flags for all reads (after some have been filtered out)
     */
    protected static JavaRDD<GATKRead> doSetPairFlags(final JavaRDD<GATKRead> reads) {

        return PSUtils.groupReadPairs(reads).flatMap(iter -> {
            boolean hasFirstMate = false;
            boolean hasSecondMate = false;
            for (final GATKRead read : iter) {
                if (read.isFirstOfPair()) hasFirstMate = true;
                else if (read.isSecondOfPair()) hasSecondMate = true;
            }
            final boolean isPaired = hasFirstMate && hasSecondMate;
            for (final GATKRead read : iter) {
                read.setIsPaired(isPaired);
            }
            return iter.iterator();
        });
    }

    @SuppressWarnings("unchecked")
    protected static JavaRDD<GATKRead> doKmerFiltering(JavaRDD<GATKRead> reads, final PipelineOptions options,
                                                       final String kmerLibPath, final byte[] mask,
                                                       final int kSize, final int countThresh) {

        reads = reads.filter(new ContainsKmerReadFilterSpark(kmerLibPath, kSize, mask, countThresh, options));
        reads.foreachPartition(read -> ContainsKmerReadFilter.closeKmerLib());
        return reads;
    }

    protected static JavaRDD<GATKRead> doBwaFilter(final JavaSparkContext ctx, final JavaRDD<GATKRead> reads,
                                                   final SAMFileHeader header, final String indexFileName,
                                                   final int minSeedLength, final int numThreads,
                                                   final int minCoverage, final int minIdentity) {

        return reads.mapPartitions(itr -> (new PSBwaFilter(indexFileName, minCoverage, minIdentity, minSeedLength, numThreads, false)).apply(itr));
    }

}
