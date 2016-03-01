package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * SparkTool to identify 63-mers in the reference that occur more than 3 times.
 */
@CommandLineProgramProperties(summary="Find the set of high copy number kmers in a reference.",
        oneLineSummary="find ref kmers with high copy number",
        programGroup = SparkProgramGroup.class)
public final class FindBadGenomicKmersSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    @VisibleForTesting static final Long MAX_KMER_FREQ = 3L;
    private static final int REF_RECORD_LEN = 10000;
    // assuming we have ~1Gb/core, we can process ~1M kmers per partition
    private static final int REF_RECORDS_PER_PARTITION = 1024*1024 / REF_RECORD_LEN;

    @Argument(doc = "file for ubiquitous kmer output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "kmer size.", fullName = "kSize", optional = true)
    private int kSize = SVConstants.KMER_SIZE;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /** Get the list of high copy number kmers in the reference, and write them to a file. */
    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final SAMFileHeader hdr = getHeaderForReads();
        SAMSequenceDictionary dict = null;
        if ( hdr != null ) dict = hdr.getSequenceDictionary();
        final PipelineOptions options = getAuthenticatedGCSOptions();
        final List<SVKmer> killList = findBadGenomicKmers(ctx, getReference(), options, dict);
        SVUtils.writeKmersFile(kSize, outputFile, options, killList);
    }

    /** Find high copy number kmers in the reference sequence */
    public static List<SVKmer> findBadGenomicKmers( final JavaSparkContext ctx,
                                                    final ReferenceMultiSource ref,
                                                    final PipelineOptions options,
                                                    final SAMSequenceDictionary readsDict ) {
        // Generate reference sequence RDD.
        final JavaRDD<byte[]> refRDD = getRefRDD(ctx, ref, options, readsDict);

        // Find the high copy number kmers
        return processRefRDD(refRDD);
    }

    /**
     * Turn a text file of overlapping records from a reference sequence into an RDD, and do a classic map/reduce:
     * Kmerize, mapping to a pair <kmer,1>, reduce by summing values by key, filter out <kmer,N> where
     * N <= MAX_KMER_FREQ, and collect the high frequency kmers back in the driver.
     */
    @VisibleForTesting static List<SVKmer> processRefRDD( final JavaRDD<byte[]> refRDD ) {
        return refRDD
                .flatMapToPair(seq ->
                        SVKmerizer.stream(seq, SVConstants.KMER_SIZE)
                                .map(kmer -> new Tuple2<>(kmer.canonical(SVConstants.KMER_SIZE),1))
                                .collect(Collectors.toCollection(() -> new ArrayList<>(seq.length))))
                .reduceByKey(Integer::sum)
                .filter(kv -> kv._2 > MAX_KMER_FREQ)
                .map(kv -> kv._1)
                .collect();
    }

    /**
     * Create an RDD from the reference sequences.
     * The reference sequences are transformed into a single, large collection of byte arrays. The collection is then
     * parallelized into an RDD.
     * Each contig that exceeds a size given by REF_RECORD_LEN is broken into a series of REF_RECORD_LEN chunks with a
     * K-1 base overlap between successive chunks. (I.e., for K=63, the last 62 bases in chunk n match the first 62
     * bases in chunk n+1) so that we don't miss any kmers due to the chunking -- we can just kmerize each record
     * independently.
     */
    private static JavaRDD<byte[]> getRefRDD( final JavaSparkContext ctx,
                                              final ReferenceMultiSource ref,
                                              final PipelineOptions options,
                                              final SAMSequenceDictionary readsDict ) {
        final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(readsDict);
        if ( dict == null ) throw new GATKException("No reference dictionary available.");

        final int effectiveRecLen = REF_RECORD_LEN - SVConstants.KMER_SIZE + 1;
        final List<byte[]> sequenceChunks = new ArrayList<>();
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final String seqName = rec.getSequenceName();
            final int seqLen = rec.getSequenceLength();
            final SimpleInterval interval = new SimpleInterval(seqName, 1, seqLen);
            try {
                final byte[] bases = ref.getReferenceBases(options, interval).getBases();
                for ( int start = 0; start < seqLen; start += effectiveRecLen ) {
                    sequenceChunks.add(Arrays.copyOfRange(bases, start, Math.min(start+REF_RECORD_LEN, seqLen)));
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't get reference sequence bases for " + interval, ioe);
            }
        }

        return ctx.parallelize(sequenceChunks, sequenceChunks.size()/REF_RECORDS_PER_PARTITION+1);
    }
}
