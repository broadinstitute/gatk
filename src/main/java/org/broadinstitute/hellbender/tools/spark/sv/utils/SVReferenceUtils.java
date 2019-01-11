package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class SVReferenceUtils {

    // TODO: 11/3/17 as pointed out in PR #3674, this method could be useful to other developers as well
    // but before that we need to have a more user friendly signature (including return type), maybe see how it is used in FindBadGenomicKmersSpark
    /**
     * Create an RDD from the reference sequences.
     * The reference sequences are transformed into a single, large collection of byte arrays.
     * The collection is then parallelized into an RDD.
     * Each contig that exceeds a size given by {@code refRecordLen} is broken into a series of {@code refRecordLen}
     * chunks with a {@code kSize} - 1 base overlap between successive chunks.
     * (I.e., for {@code kSize} = 63, the last 62 bases in chunk n match the first 62 bases in chunk n+1)
     * so that we don't miss any kmers due to the chunking -- we can just kmerize each record independently.
     */
    public static JavaRDD<byte[]> getReferenceBasesRDD(final JavaSparkContext ctx,
                                                       final int kSize,
                                                       final ReferenceMultiSparkSource ref,
                                                       final SAMSequenceDictionary dict,
                                                       final int refRecordLen,
                                                       final int refRecordsPerPartition) {
        Utils.nonNull(dict, "provided dictionary is null");
        Utils.validateArg(kSize!=0, "provided kmer size is zero");
        Utils.validateArg(refRecordLen > 0, "provided ref record length is non positive + " + refRecordLen);
        Utils.validateArg(refRecordsPerPartition > 0, "provided ref record per partition is non positive + " + refRecordsPerPartition);

        final int effectiveRecLen = refRecordLen - kSize + 1;
        final List<byte[]> sequenceChunks = new ArrayList<>();
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final String seqName = rec.getSequenceName();
            final int seqLen = rec.getSequenceLength();
            final SimpleInterval interval = new SimpleInterval(seqName, 1, seqLen);
            try {
                final byte[] bases = ref.getReferenceBases(interval).getBases();
                for ( int start = 0; start < seqLen; start += effectiveRecLen ) {
                    sequenceChunks.add(Arrays.copyOfRange(bases, start, Math.min(start+refRecordLen, seqLen)));
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't get reference sequence bases for " + interval, ioe);
            }
        }

        return ctx.parallelize(sequenceChunks, sequenceChunks.size()/refRecordsPerPartition+1);
    }

}
