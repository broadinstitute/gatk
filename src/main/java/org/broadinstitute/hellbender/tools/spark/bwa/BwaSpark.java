package org.broadinstitute.hellbender.tools.spark.bwa;

import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import com.google.common.base.Stopwatch;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.UnmodifiableIterator;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.bwa.BWANativeLibrary;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.seqdoop.hadoop_bam.BAMInputFormat;
import scala.Tuple2;

import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.function.Function;

@CommandLineProgramProperties(summary = "Runs BWA",
        oneLineSummary = "BWA on Spark",
        programGroup = SparkProgramGroup.class)
public final class BwaSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    private static final Logger log = LogManager.getLogger(BwaSpark.class);

    private static ThreadLocal<Stopwatch> stopwatchThreadLocal = new ThreadLocal<Stopwatch>() {
        @Override
        protected Stopwatch initialValue() {
            return Stopwatch.createUnstarted();
        }
    };

    private static Set<Stopwatch> stopwatches = new HashSet<>();

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Argument(doc = "the number of base pairs to send in a batch to BWA", shortName = "K",
            fullName = "fixedChunkSize", optional = true)
    private int fixedChunkSize = 100000;

    @Argument(doc = "the number of threads", shortName = "t",
            fullName = "threads", optional = true)
    private int numThreads = 1;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        Stopwatch toolStopwatch = Stopwatch.createStarted();
        // ensure reads in a pair fall in the same partition (input split), so they are processed together
        ctx.hadoopConfiguration().setBoolean(BAMInputFormat.KEEP_PAIRED_READS_TOGETHER_PROPERTY, true);

        JavaRDD<Tuple2<ShortRead, ShortRead>> shortReadPairs = getUnalignedReadPairs();

        JavaRDD<String> samLines = align(shortReadPairs);

        final SAMSequenceDictionary sequenceDictionary;
        try {
            String ref = referenceArguments.getReferenceFileName();
            FastaSequenceFile fasta = BucketUtils.isHadoopUrl(ref) ?
                    new FastaSequenceFile(IOUtils.getPath(ref), true) : new FastaSequenceFile(new File(ref), true);
            sequenceDictionary = fasta.getSequenceDictionary();
            if (sequenceDictionary == null) {
                throw new UserException("No sequence dictionary found for " + ref);
            }
        } catch (IOException e) {
            throw new GATKException("Cannot read sequence dictionary", e);
        }
        final SAMFileHeader readsHeader = new SAMFileHeader();
        readsHeader.setSortOrder(getHeaderForReads().getSortOrder());
        readsHeader.setSequenceDictionary(sequenceDictionary);

        final SAMLineParser samLineParser = new SAMLineParser(new DefaultSAMRecordFactory(), ValidationStringency.SILENT, readsHeader, null, null);
        Broadcast<SAMLineParser> samLineParserBroadcast = ctx.broadcast(samLineParser);

        JavaRDD<GATKRead> reads = samLines.map(r -> new SAMRecordToGATKReadAdapter(samLineParserBroadcast.getValue().parseLine(r)));

        try {
            ReadsSparkSink.writeReads(ctx, output, null, reads, readsHeader, shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
        } catch (IOException e) {
            throw new GATKException("Unable to write bam",e);
        }

        reportThreadTimings();
        log.info("Total time to run tool: " + toolStopwatch.elapsed(TimeUnit.SECONDS) + "s");
    }

    private JavaRDD<Tuple2<ShortRead, ShortRead>> getUnalignedReadPairs() {
        JavaRDD<GATKRead> unalignedReads = getReads();
        JavaRDD<List<GATKRead>> unalignedPairs = unalignedReads.mapPartitions(iter -> () -> Iterators.partition(iter, 2));

        return unalignedPairs.map(p -> {
            GATKRead read1 = p.get(0);
            GATKRead read2 = p.get(1);
            String name1 = read1.getName();
            String name2 = read2.getName();
            byte[] baseQualities1 = SAMUtils.phredToFastq(read1.getBaseQualities()).getBytes();
            byte[] baseQualities2 = SAMUtils.phredToFastq(read2.getBaseQualities()).getBytes();
            return new Tuple2<>(
                    new ShortRead(name1, read1.getBases(), baseQualities1),
                    new ShortRead(name2, read2.getBases(), baseQualities2));
        });
    }

    private JavaRDD<String> align(JavaRDD<Tuple2<ShortRead, ShortRead>> shortReadPairs) {
        return shortReadPairs.mapPartitions(iter -> () -> {
            BWANativeLibrary.load();

            try {
                String ref = referenceArguments.getReferenceFileName();
                File localRef = BucketUtils.isHadoopUrl(ref) ? localizeReferenceAndBwaIndexFiles(ref) : new File(ref);
                BwaIndex index = new BwaIndex(localRef);
                BwaMem mem = new BwaMem(index);
                return concat(alignChunks(mem, iter, fixedChunkSize));
            } catch (IOException e) {
                throw new GATKException("Cannot run BWA-MEM", e);
            }
        });
    }

    private File localizeReferenceAndBwaIndexFiles(String referenceFilename) throws IOException {
        Stopwatch downloadRefStopwatch = Stopwatch.createStarted();
        File localRef = File.createTempFile("ref", ".fa");
        if (!localRef.delete()) {
            throw new IOException("Cannot delete temporary file for reference: " + localRef);
        }
        Path localPath = localRef.toPath();
        Path remotePath = IOUtils.getPath(referenceFilename);
        Files.copy(remotePath, localPath);
        for (String extension : new String[] { ".amb", ".ann", ".bwt", ".pac", ".sa" }) {
            Files.copy(remotePath.resolveSibling(remotePath.getFileName() + extension),
                    localPath.resolveSibling(localPath.getFileName() + extension));
        }
        downloadRefStopwatch.stop();
        log.info("Time to download reference: " + downloadRefStopwatch.elapsed(TimeUnit.SECONDS) + "s");
        return localRef;
    }

    /**
     * Concatenates a series of {@link Iterator}s (all of the same type) into a single {@link Iterator}.
     * @param iterator an {@link Iterator} of {@link Iterator}s
     * @param <T> the type of the iterator
     * @return an {@link Iterator} over the underlying {@link Iterator}s
     */
    static <T> Iterator<T> concat(Iterator<? extends Iterable<T>> iterator) {
        return new AbstractIterator<T>() {
            Iterator<T> subIterator;
            @Override
            protected T computeNext() {
                if (subIterator != null && subIterator.hasNext()) {
                    return subIterator.next();
                }
                while (iterator.hasNext()) {
                    subIterator = iterator.next().iterator();
                    if (subIterator.hasNext()) {
                        return subIterator.next();
                    }
                }
                return endOfData();
            }
        };
    }

    /**
     * Aligns a collection ("chunk") of read pairs against the reference using BWA-MEM. Alignment is done in parallel
     * using a thread pool.
     * @param bwaMem the BWA-MEM JNI object to use to do alignment
     * @param iter the read pairs in the collection ("chunk")
     * @param fixedChunkSize the number of read bases in a chunk
     * @return an {@link Iterator} of chunks of alignment strings (in SAM format)
     */
    private Iterator<List<String>> alignChunks(final BwaMem bwaMem, Iterator<Tuple2<ShortRead, ShortRead>> iter, int fixedChunkSize) {
        return transformParallel(chunk(iter, fixedChunkSize), input -> {
            List<ShortRead> reads1 = new ArrayList<>();
            List<ShortRead> reads2 = new ArrayList<>();
            for (Tuple2<ShortRead, ShortRead> p : input) {
                reads1.add(p._1);
                reads2.add(p._2);
            }
            try {
                Stopwatch stopwatch = stopwatchThreadLocal.get();
                stopwatches.add(stopwatch);
                stopwatch.start();
                String[] alignments = bwaMem.align(reads1, reads2);
                stopwatch.stop();
                return Arrays.asList(alignments);
            } catch (IOException e) {
                throw new GATKException(e.toString());
            }
        }, numThreads);
    }

    /**
     * Breaks an {@link Iterator} over read pairs into chunks of a given size (number of bases).
     * @param iterator the read pairs
     * @param fixedChunkSize the number of read bases in a chunk
     * @return an {@link Iterator} of chunks of read pairs
     */
    private Iterator<List<Tuple2<ShortRead, ShortRead>>> chunk(final Iterator<Tuple2<ShortRead, ShortRead>> iterator, int fixedChunkSize) {
        return new UnmodifiableIterator<List<Tuple2<ShortRead, ShortRead>>>() {
            @Override
            public boolean hasNext() {
                return iterator.hasNext();
            }
            @Override
            public List<Tuple2<ShortRead, ShortRead>> next() {
                if (!hasNext()) {
                    throw new NoSuchElementException();
                }
                int size = 0;
                List<Tuple2<ShortRead, ShortRead>> list = new ArrayList<>();
                while (iterator.hasNext() && size < fixedChunkSize) {
                    Tuple2<ShortRead, ShortRead> pair = iterator.next();
                    list.add(pair);
                    size += pair._1.getBases().length;
                    size += pair._2.getBases().length;
                }
                return list;
            }
        };
    }

    /**
     * Like Guava's {@link Iterators#transform(Iterator, com.google.common.base.Function)}, but runs a fixed number
     * ({@code numThreads}) of transformations in parallel, while maintaining ordering of the output iterator.
     * This is useful if the transformations are CPU intensive.
     */
    static <F, T> Iterator<T> transformParallel(final Iterator<F> fromIterator, final Function<F, T> function, int numThreads) {
        if (numThreads == 1) { // defer to Guava for single-threaded case
            return Iterators.transform(fromIterator, new com.google.common.base.Function<F, T>() {
                @Nullable
                @Override
                public T apply(@Nullable F input) {
                    return function.apply(input);
                }
            });
        }
        // use an executor service for the multi-threaded case
        final ExecutorService executorService = Executors.newFixedThreadPool(numThreads);
        final Queue<Future<T>> futures = new LinkedList<>();
        return new AbstractIterator<T>() {
            @Override
            protected T computeNext() {
                try {
                    while (fromIterator.hasNext()) {
                        if (futures.size() == numThreads) {
                            return futures.remove().get();
                        }
                        F next = fromIterator.next();
                        Future<T> future = executorService.submit(() -> function.apply(next));
                        futures.add(future);
                    }
                    if (!futures.isEmpty()) {
                        return futures.remove().get();
                    }
                    executorService.shutdown();
                    return endOfData();
                } catch (InterruptedException | ExecutionException e) {
                    throw new GATKException("Problem running task", e);
                }
            }
        };
    }

    private static void reportThreadTimings() {
        long bwaTime = 0;
        for (Stopwatch stopwatch : stopwatches) {
            bwaTime += stopwatch.elapsed(TimeUnit.MILLISECONDS);
        }
        log.info("Time in bwa: " + (bwaTime / 1000) + "s");
    }
}
