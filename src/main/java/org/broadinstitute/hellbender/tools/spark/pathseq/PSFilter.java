package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.filters.AmbiguousBaseReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.tools.spark.pathseq.loggers.PSFilterLogger;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.ReadFilterSparkifier;
import org.broadinstitute.hellbender.tools.spark.utils.ReadTransformerSparkifier;
import org.broadinstitute.hellbender.transformers.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.illumina.IlluminaAdapterPair;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Tuple2;

import java.util.*;

/**
 * Performs PathSeq filtering steps and manages associated resources.
 */
public final class PSFilter implements AutoCloseable {

    private final JavaSparkContext ctx;
    private final PSFilterArgumentCollection filterArgs;
    private final SAMFileHeader header;

    private static final List<String> ADAPTER_SEQUENCES = CollectionUtil.makeList(
            IlluminaAdapterPair.SINGLE_END.get5PrimeAdapter(),
            IlluminaAdapterPair.SINGLE_END.get3PrimeAdapter(),
            IlluminaAdapterPair.PAIRED_END.get5PrimeAdapter(),
            IlluminaAdapterPair.PAIRED_END.get3PrimeAdapter(),
            IlluminaAdapterPair.INDEXED.get5PrimeAdapter(),
            IlluminaAdapterPair.INDEXED.get3PrimeAdapter()
    );

    private final static int REPEAT_WINDOW_SIZE_1 = 30;
    private final static int MAX_AT_CONTENT_1 = 29;
    private final static int MAX_GC_CONTENT_1 = 29;
    private final static int REPEAT_WINDOW_SIZE_2= 100;
    private final static int MAX_AT_CONTENT_2 = 87;
    private final static int MAX_GC_CONTENT_2 = 89;

    public PSFilter(final JavaSparkContext ctx, final PSFilterArgumentCollection filterArgs,
                    final SAMFileHeader header) {
        Utils.nonNull(ctx, "JavaSparkContext cannot be null");
        Utils.nonNull(filterArgs, "Filter arguments cannot be null");
        this.ctx = ctx;
        this.filterArgs = filterArgs;
        this.header = header;
    }

    @VisibleForTesting
    static JavaRDD<GATKRead> setPairFlags(final JavaRDD<GATKRead> reads, final int readsPerPartitionGuess) {
        return repartitionReadsByName(reads).mapPartitions(iter -> setPartitionUnpairedFlags(iter, readsPerPartitionGuess));
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
        if (read.isReverseStrand()) {
            SequenceUtil.reverseComplement(newRead.getBases());
            SequenceUtil.reverseQualities(newRead.getBaseQualities());
        }
        newRead.setIsUnmapped();
        newRead.setIsPaired(read.isPaired());
        if (read.isPaired()) {
            newRead.setMateIsUnmapped();
            if (read.isFirstOfPair()) {
                newRead.setIsFirstOfPair();
            } else if (read.isSecondOfPair()) {
                newRead.setIsSecondOfPair();
            }
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
     * Repartitions reads so that reads with the same name will be on the same partition
     */
    static JavaRDD<GATKRead> repartitionReadsByName(final JavaRDD<GATKRead> reads) {
        return repartitionReadsByName(reads, reads.getNumPartitions());
    }

    static JavaRDD<GATKRead> repartitionReadsByName(final JavaRDD<GATKRead> reads, final int numPartitions) {
        //Shuffle reads into partitions by read name hash code
        return reads.mapToPair(read -> new Tuple2<>(read.getName(), read))
                .partitionBy(new HashPartitioner(numPartitions))
                .map(Tuple2::_2);
    }

    /**
     * Maps partition to a Tuple of two Lists, the first containing the paired reads, the second containing unpaired
     */
    static Tuple2<List<GATKRead>, List<GATKRead>> getPairedAndUnpairedLists(final Iterator<GATKRead> iter, final int readsPerPartitionGuess) {
        //Find the paired and unpaired reads by scanning the partition for repeated names
        final ArrayList<GATKRead> pairedReadsList = new ArrayList<>(readsPerPartitionGuess);
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
        pairedReadsList.trimToSize();

        return new Tuple2<>(pairedReadsList, unpairedReadsList);
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
        final long hashForward = SVUtils.fnvByteArray64(bases);
        SequenceUtil.reverseComplement(bases);
        final long hashReverse = SVUtils.fnvByteArray64(bases);
        return new Tuple2<>(Math.min(hashForward, hashReverse), read);
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
                                         final int minIdentity) {

        return reads.mapPartitions(itr -> (new PSBwaFilter(indexFileName, minIdentity, minSeedLength, numThreads, false)).apply(itr));
    }

    /**
     * Main PathSeq filtering method. See PathSeqFilterSpark for an overview.
     * Returns a tuple containing the paired reads and unpaired reads as separate RDDs.
     * If metricsFile is null, read count metrics will not be collected.
     */
    public Tuple2<JavaRDD<GATKRead>, JavaRDD<GATKRead>> doFilter(JavaRDD<GATKRead> reads, final PSFilterLogger filterLogger) {

        Utils.nonNull(reads, "Input reads cannot be null");
        reads = PSUtils.primaryReads(reads);
        filterLogger.logPrimaryReads(reads);

        if (filterArgs.alignedInput) {
            reads = reads.filter(new ReadFilterSparkifier(new HostAlignmentReadFilter(filterArgs.minIdentity)));
        }
        filterLogger.logReadsAfterPrealignedHostFilter(reads);

        //Clear alignment data from the reads
        reads = clearAllAlignments(reads, header);

        //Remove /1 and /2 from read names
        reads = reads.map(new ReadTransformerSparkifier(new StripMateNumberTransformer()));

        if (!filterArgs.skipFilters) {

            //Adapter trimming
            reads = reads.map(new ReadTransformerSparkifier(new AdapterTrimTransformer(filterArgs.maxAdapterMismatches, filterArgs.minAdapterLength, ADAPTER_SEQUENCES)));

            //Apply simple repeat masking
            //See "Low-complexity DNA and simple repeats" at http://www.repeatmasker.org/webrepeatmaskerhelp.html
            reads = reads.map(new ReadTransformerSparkifier(new SimpleRepeatMaskTransformer(MAX_AT_CONTENT_1, MAX_GC_CONTENT_1, REPEAT_WINDOW_SIZE_1)));
            reads = reads.map(new ReadTransformerSparkifier(new SimpleRepeatMaskTransformer(MAX_AT_CONTENT_2, MAX_GC_CONTENT_2, REPEAT_WINDOW_SIZE_2)));

            //Apply DUST masking
            reads = reads.map(new ReadTransformerSparkifier(new DUSTReadTransformer(filterArgs.dustMask, filterArgs.dustW, filterArgs.dustT)));

            //Apply base quality hard clipping
            reads = reads.map(new ReadTransformerSparkifier(new BaseQualityClipReadTransformer(filterArgs.readTrimThresh)));

            //Filter reads with less than minReadLength bases
            reads = reads.filter(new ReadFilterSparkifier(new ReadLengthReadFilter(filterArgs.minReadLength, Integer.MAX_VALUE)));

            //Change low-quality bases to 'N'
            reads = reads.map(new ReadTransformerSparkifier(new BaseQualityReadTransformer(filterArgs.qualPhredThresh)));

            //Filter reads with too many 'N's
            reads = reads.filter(new ReadFilterSparkifier(new AmbiguousBaseReadFilter(filterArgs.maxAmbiguousBases)));
        }
        filterLogger.logReadsAfterQualityFilter(reads);

        //Kmer filtering
        if (filterArgs.kmerFilePath != null) {
            reads = doKmerFiltering(reads, filterArgs.kmerFilePath, filterArgs.hostKmerThresh);
        }

        //Redistribute reads
        if (!filterArgs.skipPreBwaRepartition) {
            reads = repartitionReadsByName(reads);
        }

        //Bwa host alignment filtering
        if (filterArgs.indexImageFile != null) {
            reads = doBwaFilter(reads, filterArgs.indexImageFile, filterArgs.minSeedLength,
                    filterArgs.bwaThreads, filterArgs.minIdentity);
        }
        filterLogger.logReadsAfterHostFilter(reads);

        //Filter duplicates
        if (filterArgs.filterDuplicates) {
            reads = setPairFlags(reads, filterArgs.filterReadsPerPartition);
            reads = filterDuplicateSequences(reads);
        }
        filterLogger.logReadsAfterDeduplication(reads);

        //Sets pairedness flags properly
        reads = setPairFlags(reads, filterArgs.filterReadsPerPartition);
        reads = clearAllAlignments(reads, header);

        //Unset paired read flags for reads that are not paired
        final PSPairedUnpairedSplitterSpark splitter = new PSPairedUnpairedSplitterSpark(reads, filterArgs.filterReadsPerPartition, false);
        final JavaRDD<GATKRead> pairedReads = splitter.getPairedReads();
        final JavaRDD<GATKRead> unpairedReads = splitter.getUnpairedReads();
        filterLogger.logFinalPairedReads(pairedReads);

        return new Tuple2<>(pairedReads, unpairedReads);
    }

    /**
     * After doFilter(), this should be run after a Spark action (e.g. write bam) has been invoked on both output RDDs
     */
    public void close() {
        BwaMemIndexCache.closeAllDistributedInstances(ctx);
        ContainsKmerReadFilterSpark.closeAllDistributedInstances(ctx);
    }

}
