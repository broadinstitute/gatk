package org.broadinstitute.hellbender.tools.spark.bwa;

import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Iterators;
import com.google.common.collect.UnmodifiableIterator;
import htsjdk.samtools.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BWANativeLibrary;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.seqdoop.hadoop_bam.BAMInputFormat;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * A reusable core functionality of running BWA on Spark.
 * The principal way to interact with the engine is like this:
 <code>
 BwaSparkEngine engine = new BwaSparkEngine(numThreads, fixedChunkSize, referenceFileName);
 SAMFileHeader readsHeader = engine.makeSamFileHeader(getHeaderForReads().getSortOrder());
 JavaRDD<GATKRead> reads = engine.alignWithBWA(ctx, unalignedReads, readsHeader);
 </code>
 */
public final class BwaSparkEngine implements Serializable{

    private static final Logger log = LogManager.getLogger(BwaSparkEngine.class);

    private static final long serialVersionUID = 1L;
    private final int numThreads;
    private final int fixedChunkSize;
    private final String referenceFileName;

    /**
     * @param numThreads number of threads to use (at least 1)
     * @param fixedChunkSize - the number of base pairs to send in a batch to BWA (at least 1) - corresponds to -K in BWA
     * @param referenceFileName
     */
    public BwaSparkEngine(final int numThreads, final int fixedChunkSize, final String referenceFileName){
        Utils.validateArg(numThreads >= 1, "num threads");
        Utils.validateArg(fixedChunkSize >= 1, "fixedChunkSize");
        Utils.nonNull(referenceFileName, "referenceFileName");
        this.numThreads = numThreads;
        this.fixedChunkSize = fixedChunkSize;
        this.referenceFileName = referenceFileName;
    }

    public JavaRDD<GATKRead> alignWithBWA(final JavaSparkContext ctx, final JavaRDD<GATKRead> unalignedReads, final SAMFileHeader readsHeader) {
        //Note: SparkContext is not serializable so we don't store it in the engine and set this property here. Setting it multiple times is fine.
        // ensure reads in a pair fall in the same partition (input split), so they are processed together
        ctx.hadoopConfiguration().setBoolean(BAMInputFormat.KEEP_PAIRED_READS_TOGETHER_PROPERTY, true);

        final JavaRDD<Tuple2<ShortRead, ShortRead>> shortReadPairs = convertToUnalignedReadPairs(unalignedReads);
        final JavaRDD<String> samLines = align(shortReadPairs);
        final Broadcast<SAMFileHeader>  readsHeaderBroadcast = ctx.broadcast(readsHeader);
        return samLines.mapPartitions(lineIterator -> {
            //Note: The parser is stateful and not thread safe.
            // So we reuse the parser over the whole partition and call it from one thread only.
            final SAMLineParser samLineParser = new SAMLineParser(new DefaultSAMRecordFactory(), ValidationStringency.SILENT, readsHeaderBroadcast.getValue(), null, null);
            final List<GATKRead> reads = new ArrayList<>();
            while (lineIterator.hasNext()) {
                final String samLine = lineIterator.next();
                reads.add(new SAMRecordToGATKReadAdapter(samLineParser.parseLine(samLine)));
            }
            return reads.iterator();
        });
    }

    private JavaRDD<Tuple2<ShortRead, ShortRead>> convertToUnalignedReadPairs(final JavaRDD<GATKRead> unalignedReads) {
        final JavaRDD<List<GATKRead>> unalignedPairs = unalignedReads.mapPartitions(iter -> Iterators.partition(iter, 2));

        return unalignedPairs.map(p -> {
            final GATKRead read1 = p.get(0);
            final GATKRead read2 = p.get(1);
            final String name1 = read1.getName();
            final String name2 = read2.getName();
            final byte[] baseQualities1 = SAMUtils.phredToFastq(read1.getBaseQualities()).getBytes();
            final byte[] baseQualities2 = SAMUtils.phredToFastq(read2.getBaseQualities()).getBytes();
            return new Tuple2<>(
                    new ShortRead(name1, read1.getBases(), baseQualities1),
                    new ShortRead(name2, read2.getBases(), baseQualities2));
        });
    }

    private JavaRDD<String> align(final JavaRDD<Tuple2<ShortRead, ShortRead>> shortReadPairs) {
        return shortReadPairs.mapPartitions(iter -> {
            BWANativeLibrary.load();

            try {
                final File localRef = BucketUtils.isHadoopUrl(referenceFileName) ? localizeReferenceAndBwaIndexFiles(referenceFileName) : new File(referenceFileName);
                final BwaIndex index = new BwaIndex(localRef);
                final BwaMem mem = new BwaMem(index);
                return Utils.concatIterators(alignChunks(mem, iter));
            } catch (final IOException e) {
                throw new GATKException("Cannot run BWA-MEM", e);
            }
        });
    }


    public static File localizeReferenceAndBwaIndexFiles(final String referenceFilename) throws IOException {
        final Stopwatch downloadRefStopwatch = Stopwatch.createStarted();
        final File localRef = File.createTempFile("referenceFileName", ".fa");
        if (!localRef.delete()) {
            throw new IOException("Cannot delete temporary file for reference: " + localRef);
        }
        final Path localPath = localRef.toPath();
        final Path remotePath = IOUtils.getPath(referenceFilename);
        Files.copy(remotePath, localPath);
        for (final String extension : new String[] { ".amb", ".ann", ".bwt", ".pac", ".sa" }) {
            Files.copy(remotePath.resolveSibling(remotePath.getFileName() + extension),
                    localPath.resolveSibling(localPath.getFileName() + extension));
        }
        downloadRefStopwatch.stop();
        log.info("Time to download reference: " + downloadRefStopwatch.elapsed(TimeUnit.SECONDS) + "s");
        return localRef;
    }

    /**
     * Aligns a collection ("chunk") of read pairs against the reference using BWA-MEM. Alignment is done in parallel
     * using a thread pool.
     * @param bwaMem the BWA-MEM JNI object to use to do alignment
     * @param iter the read pairs in the collection ("chunk")
     * @return an {@link Iterator} of chunks of alignment strings (in SAM format)
     */
    private Iterator<List<String>> alignChunks(final BwaMem bwaMem, final Iterator<Tuple2<ShortRead, ShortRead>> iter) {
        return Utils.transformParallel(chunk(iter), input -> {
            final List<ShortRead> reads1 = new ArrayList<>();
            final List<ShortRead> reads2 = new ArrayList<>();
            for (final Tuple2<ShortRead, ShortRead> p : input) {
                reads1.add(p._1);
                reads2.add(p._2);
            }
            try {
                final String[] alignments = bwaMem.align(reads1, reads2);
                return Arrays.asList(alignments);
            } catch (final IOException e) {
                throw new GATKException(e.toString());
            }
        }, numThreads);
    }


    /**
     * Breaks an {@link Iterator} over read pairs into chunks of a given size (number of bases).
     * @param iterator the read pairs
     * @return an {@link Iterator} of chunks of read pairs
     */
    private Iterator<List<Tuple2<ShortRead, ShortRead>>> chunk(final Iterator<Tuple2<ShortRead, ShortRead>> iterator) {
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
                final List<Tuple2<ShortRead, ShortRead>> list = new ArrayList<>();
                while (iterator.hasNext() && size < fixedChunkSize) {
                    final Tuple2<ShortRead, ShortRead> pair = iterator.next();
                    list.add(pair);
                    size += pair._1.getBases().length;
                    size += pair._2.getBases().length;
                }
                return list;
            }
        };
    }

    /**
     * Creates a header for the SAM/BAM file that will be written.
     */
    public SAMFileHeader makeHeaderForOutput(final SAMFileHeader inputHeader, final SAMSequenceDictionary refDictionary) {
        Utils.nonNull(inputHeader);
        if (inputHeader.getSequenceDictionary() != null && !inputHeader.getSequenceDictionary().isEmpty()){
            return inputHeader;
        }
        if (refDictionary == null) {
            throw new UserException("No sequence dictionary found in the reference " + referenceFileName);
        }
        final SAMFileHeader readsHeader = inputHeader.clone();//make a copy and set the dictionary from the reference.
        readsHeader.setSequenceDictionary(refDictionary);
        return readsHeader;
    }
}
