package org.broadinstitute.hellbender.utils.nio;

import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * NioBam holds a reference to a BAM file on Google Cloud Storage, and can give you
 * an RDD with the reads inside of it.
 *
 * The constructor will open the file to make sure it exists and we have access to it.
 * This is preferable to waiting until we're running on the cloud.
 *
 * Although we don't expect you to move these objects across computers, you can.
 */
public class NioBam implements Serializable {
    protected static final Logger logger = LogManager.getLogger(NioBam.class);

    private static final long serialVersionUID = 1L;
    private final String bam;
    private final String index;
    private transient byte[] indexCache;

    /** Checks the files exists, then stores them. **/
    public NioBam(String gcsFilename, String indexGcsFilename) {
        try {
            this.bam = gcsFilename;
            this.index = indexGcsFilename;
            init();
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile("Could not read file " + gcsFilename, e);
        }
    }

    /** Finds the index file, then calls NioBam(bam, index). **/
    public NioBam(String gcsFilename) {
        try {
            String indexFilename = gcsFilename + ".bai";
            if ( !Files.exists(IOUtils.getPath(indexFilename)) ) {
                int i = gcsFilename.lastIndexOf('.');
                if ( i >= 0 ) {
                    indexFilename = gcsFilename.substring(0, i) + ".bai";
                }
            }
            this.bam = gcsFilename;
            this.index = indexFilename;
            init();
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile("Could not read file " + gcsFilename, e);
        }
    }

    private void init() throws FileNotFoundException {
        Path bamPath = IOUtils.getPath(bam);
        Path bamIndexPath = IOUtils.getPath(index);
        if (!Files.exists(bamPath)) {
            throw new FileNotFoundException(bamPath.toString());
        }
        if (!Files.exists(bamIndexPath)) {
            throw new FileNotFoundException(bamIndexPath.toString());
        }
    }

    /** Parses the BAM file into SAMRecords. Will be distributed onto at least 'numPartitions' partitions. **/
    public JavaRDD<SAMRecord> getReads(JavaSparkContext ctx, int numPartitions) {
        try {
            Path bamPath = IOUtils.getPath(bam);
            ChannelAsSeekableStream bamOverNIO = new ChannelAsSeekableStream(Files.newByteChannel(bamPath), bamPath.toString());
            final byte[] index = getIndex();
            SeekableStream indexInMemory = new ByteArraySeekableStream(index);

            SamReader bam3 = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.LENIENT)
                    .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
                    .open(SamInputResource.of(bamOverNIO).index(indexInMemory));
            List<QueryInterval[]> chunks = getAllChunksBalanced(bam3, numPartitions, numPartitions);

            // Ideally we'd get exactly the number of chunks the user is asking for, but until then...
            logger.info("We got: " + chunks.size() + " batches of chunks.");

            return ctx.parallelize(chunks, chunks.size()).flatMap(qis -> new ReadsIterable(bam, index, qis).iterator());
        }
        catch ( IOException e ) {
            throw new GATKException("I/O error loading reads", e);
        }
    }

    private synchronized byte[] getIndex() throws IOException {
        if (null!= indexCache) {
            return indexCache;
        }
        indexCache = Files.readAllBytes(IOUtils.getPath(index));
        return indexCache;
    }

    private static class SizedQueryInterval {
        public QueryInterval queryInterval;
        public long size;

        public SizedQueryInterval(QueryInterval queryInterval, long size) {
            this.queryInterval = queryInterval;
            this.size = size;
        }
    }

    // this isn't very good yet, ideally we want just this number of query intervals, not per-contig.
    private static List<QueryInterval[]> getAllChunksBalanced(SamReader bam, int countPerContig, int totalCount) {
        List<SizedQueryInterval> ret = new ArrayList<>();
        SAMFileHeader header = bam.getFileHeader();
        for (SAMSequenceRecord s : header.getSequenceDictionary().getSequences()) {
            ret.addAll(getChunksBalanced(bam, s.getSequenceIndex(), countPerContig));
        }
        return rebalanceChunks(ret, totalCount);
    }

    // Given a list of SizedQueryInterval, put them into "retCount" batches of approximately the
    // same size.
    // Try to keep neighboring chunks together, and keep the computation fast (so it's definitely
    // not returning anything near the optimal solution).
    private static List<QueryInterval[]> rebalanceChunks(List<SizedQueryInterval> chunks, int retCount) {
        long shortestBatch = 999999999999999999L;
        long longestBatch = 0L;
        long totalSoFar = 0;
        long batchesSoFar = 0;
        List<QueryInterval[]> ret = new ArrayList<>();
        long totalSize = 0;
        for (SizedQueryInterval c : chunks) {
            totalSize += c.size;
        }
        long tgtSize = totalSize / retCount;
        // first off, we grab all of the chunks that are larger than target size
        List<SizedQueryInterval> smallChunks = new ArrayList<>();
        for (SizedQueryInterval c : chunks) {
            if (c.size >= tgtSize) {
                if (c.size > longestBatch) longestBatch = c.size;
                if (c.size < shortestBatch) shortestBatch = c.size;
                List<QueryInterval> batch = new ArrayList<>(1);
                batch.add(c.queryInterval);
                ret.add(optimize(batch));
                totalSoFar += c.size;
                batchesSoFar++;
            } else {
                smallChunks.add(c);
            }
        }
        // next, add up the little chunks until they are large enough
        long sizeSoFar = 0;
        List<QueryInterval> batch = new ArrayList<>();
        for (SizedQueryInterval c : smallChunks) {
            long newSize = sizeSoFar + c.size;
            // adjust the target batch size based on how many we have left.
            long adjustedTargetSize = tgtSize;
            if (batchesSoFar < retCount) {
                adjustedTargetSize = (totalSize - totalSoFar) / (retCount - batchesSoFar);
            }
            if (newSize >= adjustedTargetSize) {
                if ((sizeSoFar + newSize) / 2.0 <= adjustedTargetSize) {
                    // We're closer to the target by including the latest guy
                    if (newSize > longestBatch) longestBatch = newSize;
                    if (newSize < shortestBatch) shortestBatch = newSize;
                    batch.add(c.queryInterval);
                    ret.add(optimize(batch));
                    batch = new ArrayList<>();
                    totalSoFar += newSize;
                    batchesSoFar++;
                    sizeSoFar = 0;
                } else {
                    // We're closer to the target by excluding the latest guy
                    if (sizeSoFar > longestBatch) longestBatch = sizeSoFar;
                    if (sizeSoFar < shortestBatch) shortestBatch = sizeSoFar;
                    ret.add(optimize(batch));
                    batch = new ArrayList<>();
                    batch.add(c.queryInterval);
                    totalSoFar += sizeSoFar;
                    batchesSoFar++;
                    sizeSoFar = c.size;
                }
            } else {
                batch.add(c.queryInterval);
                sizeSoFar = newSize;
            }
        }
        if (!batch.isEmpty()) {
            ret.add(optimize(batch));
        }
        logger.info("Batch size range: " + shortestBatch + " to " + longestBatch);
        return ret;
    }

    private static QueryInterval[] optimize(List<QueryInterval> list) {
        QueryInterval[] array = new QueryInterval[0];
        return QueryInterval.optimizeIntervals(list.toArray(array));
    }

    private static List<SizedQueryInterval> getChunksBalanced(SamReader bam, int sequenceIndex, int retCount) {
        List<SizedQueryInterval> ret = new ArrayList<>();
        BAMIndex index = bam.indexing().getIndex();
        SAMFileHeader header = bam.getFileHeader();
        SAMSequenceRecord s = header.getSequence(sequenceIndex);
        long totalLength = chunksLength(getChunks(index, sequenceIndex, 1, s.getSequenceLength() + 1));
        if (totalLength == 0) {
            return ret;
        }
        int sofar = 0;
        long targetLength = totalLength / retCount;
        int end = s.getSequenceLength();
        int step = s.getSequenceLength() / (100 * retCount);
        if (step < 1) step = 1;
        int start = 1;
        for (int j = step; j < end; j += step) {
            if (j > end) j = end;
            List<Chunk> candidate = getChunks(index, sequenceIndex, start, j);
            long size = chunksLength(candidate);
            if (size < targetLength) {
                // not big enough yet
                continue;
            }
            if (size > targetLength * 2) {
                // too large, search for a good separation point
                // TODO
            }
            // good, emit.
            ret.add(new SizedQueryInterval(new QueryInterval(sequenceIndex, start, j + 1), size));
            start = j;
            sofar += size;
            if (ret.size() < retCount) {
                targetLength = (totalLength - sofar) / (retCount - ret.size());
            } else {
                targetLength = totalLength / retCount;
            }

        }
        return ret;
    }


    private static List<Chunk> getChunks(BAMIndex index, int sequenceIndex, int start, int endExcluded) {
        if (endExcluded <= start) return new ArrayList<>();
        BAMFileSpan span = index.getSpanOverlapping(sequenceIndex, start, endExcluded - 1);
        if (null == span) return new ArrayList<>();
        return span.getChunks();
    }

    private static long chunksLength(List<Chunk> chunks) {
        long totalLength = 0;
        for (Chunk c : chunks) {
            totalLength += chunkSize(c);
        }
        return totalLength;
    }

    private static long chunkSize(Chunk c) {
        long start = c.getChunkStart() >> 16;
        long end = (c.getChunkEnd() >> 16) + (c.getChunkEnd() & 0xffff);
        return (end - start);
    }

}
