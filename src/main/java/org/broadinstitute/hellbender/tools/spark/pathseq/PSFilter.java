package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.filters.AmbiguousBaseReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.ReadFilterSparkifier;
import org.broadinstitute.hellbender.tools.spark.utils.ReadTransformerSparkifier;
import org.broadinstitute.hellbender.transformers.BaseQualityClipReadTransformer;
import org.broadinstitute.hellbender.transformers.BaseQualityReadTransformer;
import org.broadinstitute.hellbender.transformers.DUSTReadTransformer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexSingleton;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Tuple2;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 * Performs PathSeq filtering steps and manages associated resources.
 */
public final class PSFilter implements AutoCloseable {

    private final JavaSparkContext ctx;
    private final PSFilterArgumentCollection filterArgs;
    private final MetricsState metricsState;
    private JavaRDD<GATKRead> reads;
    private final SAMFileHeader header;

    public PSFilter(final JavaSparkContext ctx, final PSFilterArgumentCollection filterArgs, final JavaRDD<GATKRead> inputReads,
                    final SAMFileHeader header) {
        Utils.nonNull(ctx);
        Utils.nonNull(filterArgs);
        Utils.nonNull(inputReads);
        this.ctx = ctx;
        this.filterArgs = filterArgs;
        this.metricsState = initializeMetics(filterArgs.metricsFileUri);
        this.reads = inputReads;
        this.header = header;
    }

    @VisibleForTesting
    static JavaRDD<GATKRead> setPairFlags(final JavaRDD<GATKRead> reads, final int readsPerPartitionGuess) {
        return repartitionPairedReads(reads).mapPartitions(iter -> setPartitionUnpairedFlags(iter, readsPerPartitionGuess));
    }

    private static JavaRDD<GATKRead> clearAllAlignments(final JavaRDD<GATKRead> reads, final SAMFileHeader header) {
        return reads.map(read -> clearReadAlignment(read, header));
    }

    /**
     * Returns input read with alignment-related info cleared
     */
    private static GATKRead clearReadAlignment(final GATKRead read, final SAMFileHeader header) {
        final GATKRead newRead = new SAMRecordToGATKReadAdapter(new SAMRecord(header));
        newRead.setName(read.getName());
        newRead.setBases(read.getBases());
        newRead.setBaseQualities(read.getBaseQualities());
        newRead.setIsUnmapped();
        newRead.setIsPaired(read.isPaired());
        if (read.isFirstOfPair()) {
            newRead.setIsFirstOfPair();
        } else if (read.isSecondOfPair()) {
            newRead.setIsSecondOfPair();
        }
        final String readGroup = read.getReadGroup();
        if (readGroup != null) {
            newRead.setAttribute(SAMTag.RG.name(), readGroup);
        }
        return newRead;
    }

    /**
     * Sets proper pairedness flags
     */
    private static Iterator<GATKRead> setPartitionUnpairedFlags(final Iterator<GATKRead> iter, final int readsPerPartitionGuess) {
        final Tuple2<List<GATKRead>, List<GATKRead>> lists = getPairedAndUnpairedLists(iter, readsPerPartitionGuess);
        final List<GATKRead> pairedReadsList = lists._1;
        final List<GATKRead> unpairedReadsList = lists._2;
        for (final GATKRead unpairedRead : unpairedReadsList) {
            unpairedRead.setIsPaired(false);
        }
        final List<GATKRead> newPartitionList = new ArrayList<>(pairedReadsList.size() + unpairedReadsList.size());
        newPartitionList.addAll(pairedReadsList);
        newPartitionList.addAll(unpairedReadsList);
        return newPartitionList.iterator();
    }

    /**
     * Maps partition to a Tuple of two Lists, the first containing the paired reads, the second containing unpaired
     */
    static JavaRDD<GATKRead> repartitionPairedReads(final JavaRDD<GATKRead> reads) {
        //Shuffle reads into partitions by read name hash code
        return reads.mapToPair(read -> new Tuple2<>(read.getName(), read))
                .partitionBy(new HashPartitioner(reads.getNumPartitions()))
                .map(Tuple2::_2);
    }

    /**
     * Maps partition to a Tuple of two Lists, the first containing the paired reads, the second containing unpaired
     */
    static Tuple2<List<GATKRead>, List<GATKRead>> getPairedAndUnpairedLists(final Iterator<GATKRead> iter, final int readsPerPartitionGuess) {
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

        return new Tuple2<>(pairedReadsListResized, unpairedReadsList);
    }

    /**
     * Preferentially filters unpaired reads, when possible. Assumes reads have pairedness flags set properly
     */
    @VisibleForTesting
    static JavaRDD<GATKRead> filterDuplicateSequences(final JavaRDD<GATKRead> reads) {
        return reads.mapToPair(PSFilter::canonicalizeRead)
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
    @VisibleForTesting
    static Tuple2<Long, GATKRead> canonicalizeRead(final GATKRead read) {
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
            for (int j = 0; j < 8 && i + j < toHash.length; j++) {
                val = (val << 8) | toHash[i + j];
            }
            start = SVUtils.fnvLong64(start, val);
        }
        return start;
    }

    @SuppressWarnings("unchecked")
    private static JavaRDD<GATKRead> doKmerFiltering(final JavaRDD<GATKRead> reads, final String kmerLibPath,
                                                       final int countThresh) {

        return reads.filter(new ContainsKmerReadFilterSpark(kmerLibPath, countThresh));
    }

    @VisibleForTesting
    static JavaRDD<GATKRead> doBwaFilter(final JavaRDD<GATKRead> reads,
                                                   final String indexFileName,
                                                   final int minSeedLength, final int numThreads,
                                                   final int minCoverage, final int minIdentity) {

        return reads.mapPartitions(itr -> (new PSBwaFilter(indexFileName, minCoverage, minIdentity, minSeedLength, numThreads, false)).apply(itr));
    }

    /**
     * Main PathSeq filtering method. See PathSeqFilterSpark for an overview.
     * Returns a tuple containing the paired reads and unpaired reads as separate RDDs.
     */
    public Tuple2<JavaRDD<GATKRead>, JavaRDD<GATKRead>> doFilter() {

        recordReadCountMetric(reads, "input");
        reads = PSUtils.primaryReads(reads);
        recordReadCountMetric(reads, "primary_reads");

        if (filterArgs.alignedInput) {
            reads = reads.filter(new ReadFilterSparkifier(new HostAlignmentReadFilter(filterArgs.minCoverage, filterArgs.minIdentity)));
            recordReadCountMetric(reads, "prealigned_filter");
        }

        //Clear alignment data from the reads
        reads = clearAllAlignments(reads, header);

        if (!filterArgs.skipFilters) {

            //Apply DUST masking
            reads = reads.map(new ReadTransformerSparkifier(new DUSTReadTransformer(filterArgs.dustMask, filterArgs.dustW, filterArgs.dustT)));

            //Apply base quality hard clipping
            reads = reads.map(new ReadTransformerSparkifier(new BaseQualityClipReadTransformer(filterArgs.readTrimThresh)));

            //Filter reads with less than minReadLength bases
            reads = reads.filter(new ReadFilterSparkifier(new ReadLengthReadFilter(filterArgs.minReadLength, Integer.MAX_VALUE)));
            recordReadCountMetric(reads, "quality_1");

            //Change low-quality bases to 'N'
            reads = reads.map(new ReadTransformerSparkifier(new BaseQualityReadTransformer(filterArgs.qualPhredThresh)));

            //Filter reads with too many 'N's
            reads = reads.filter(new ReadFilterSparkifier(new AmbiguousBaseReadFilter(filterArgs.fracNThreshold)));
            recordReadCountMetric(reads, "quality_2");
        }

        //Kmer filtering
        if (filterArgs.kmerLibPath != null) {
            reads = doKmerFiltering(reads, filterArgs.kmerLibPath, filterArgs.hostKmerThresh);
            recordReadCountMetric(reads, "kmer");
        }

        //Bwa host alignment filtering
        if (filterArgs.indexImageFile != null) {
            reads = doBwaFilter(reads, filterArgs.indexImageFile, filterArgs.minSeedLength,
                    filterArgs.bwaThreads, filterArgs.minCoverage, filterArgs.minIdentity);
            recordReadCountMetric(reads, "bwa");
        }

        //Filter duplicates
        if (filterArgs.filterDuplicates) {
            reads = setPairFlags(reads, filterArgs.readsPerPartition);
            reads = filterDuplicateSequences(reads);
            recordReadCountMetric(reads, "duplicates");
        }

        //Unset paired read flags for reads that are not paired
        final PSPairedUnpairedSplitterSpark splitter = new PSPairedUnpairedSplitterSpark(reads, filterArgs.readsPerPartition);
        final JavaRDD<GATKRead> pairedReads = splitter.getPairedReads();
        final JavaRDD<GATKRead> unpairedReads = splitter.getUnpairedReads();
        recordReadCountMetric(reads, "split_paired");

        return new Tuple2<>(pairedReads, unpairedReads);
    }

    /**
     * Open metrics file, to which the number of reads after each stage will be written.
     */
    private MetricsState initializeMetics(final String metricsFileUri) {
        if (metricsFileUri == null) return null;
        try {
            return new MetricsState(metricsFileUri);
        } catch (IOException e) {
            throw new GATKException("Could not open metrics file " + metricsFileUri, e);
        }
    }

    /**
     * Records number of reads in the metrics file. Note this invokes a count action on the RDD and therefore
     * may substantially reduce performance.
     */
    private void recordReadCountMetric(final JavaRDD<GATKRead> reads, final String name) {
        if (metricsState == null) return;
        try {
            metricsState.writeReadCount(reads, name);
        } catch (final IOException e) {
            throw new GATKException("Could not write to metrics file", e);
        }
    }

    /**
     * After doFilter(), this should be run after a Spark action (e.g. write bam) has been invoked on both output RDDs
     */
    public void close() {
        reads.foreachPartition(read -> ContainsKmerReadFilter.closeKmerLib());
        BwaMemIndexSingleton.closeAllDistributedInstances(ctx);
        if (metricsState != null) {
            try {
                metricsState.close();
            } catch (final IOException e) {
                throw new GATKException("Could not close metrics output stream");
            }
        }
    }

    /**
     * Tracks elapsed time and output stream for the metrics file
     */
    private static final class MetricsState implements AutoCloseable {
        public final long initialTimeMillis;
        public final OutputStream metricsOutputStream;

        public MetricsState(final String metricsFileUri) throws IOException {
            initialTimeMillis = System.currentTimeMillis();
            metricsOutputStream = BucketUtils.createFile(metricsFileUri);
            final String outputString = "step\treads\tseconds\n";
            metricsOutputStream.write(outputString.getBytes(Charset.defaultCharset()));
        }

        public void writeReadCount(JavaRDD<GATKRead> reads, final String name) throws IOException {
            final long currentTimeMillis = System.currentTimeMillis();
            final double elapsedTimeSeconds = (currentTimeMillis - initialTimeMillis) / 1000.0;
            NumberFormat formatter = new DecimalFormat("#0.00");
            final String outputString = name + "\t" + reads.count() + "\t" + formatter.format(elapsedTimeSeconds) + "\n";
            metricsOutputStream.write(outputString.getBytes(Charset.defaultCharset()));
            metricsOutputStream.flush();
        }

        public void close() throws IOException {
            metricsOutputStream.close();
        }
    }
}
