package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.runners.DataflowPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.base.Stopwatch;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.dataflow.DoFnWLog;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.File;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * see ShardedReadsDatasource.get
 */
public final class ShardedReadsDatasource {
    private final static Logger logger = LogManager.getLogger(ShardedReadsDatasource.class);

    /**
     * Create a {@link PCollection< KV <String,Iterable<GATKRead>>>} containing all the reads starting in the given intervals and their mates,
     * grouped by the shard they start in.
     * The shards are strings as per SimpleInterval's ToString syntax.
     * Reads that are malformed or unmapped are ignored.
     *
     * Note: for performance reasons, this will give you SAMRecords that have no header. You have to put it back in.
     *
     * @param bam the file to open, on GCS or local disk.
     * @param intervals the subset of the file you want returned.
     * @param readShardSize shards of this size are assigned to individual machines (streamed in from GCS; larger is better).
     * @param  outputShardSize reads over this many bases are output at a time (held all together in RAM; too large and you'll go boom).
     */
    public static PCollection<KV<String,Iterable<GATKRead>>> get(String bam, List<SimpleInterval> intervals, int readShardSize, int outputShardSize, Pipeline pipeline, GenomicsFactory.OfflineAuth auth) {
        final boolean cloudStorageUrl = BucketUtils.isCloudStorageUrl(bam);
        final boolean hadoopUrl = BucketUtils.isHadoopUrl(bam);

        List<SimpleInterval> shards = cutToShards(intervals, readShardSize);
        // so if there are empty ranges, they are spread evenly across the workers.
        Collections.shuffle(shards);
        logger.info("Input range sharded into "+shards.size()+" intervals.");
        List<KV<String,Void>> intervalKeys;
        intervalKeys = shards.stream().map(s-> KV.<String,Void>of(s.toString(), null)).collect(Collectors.<KV<String,Void>>toList());

        return pipeline
            // ship shards to workers
            .apply(Create.<SimpleInterval>of(shards)).setName("shards")
                // for each shard, load the corresponding reads
            .apply(ParDo.named("getGroupedReadPCollection").of(new DoFnWLog<SimpleInterval, KV<String, Iterable<GATKRead>>>("getGroupedReadPCollection") {
                private static final long serialVersionUID = 1L;
                private transient Storage.Objects storageClient;
                private transient Stopwatch querying;
                private transient Stopwatch loading;
                private transient Stopwatch choosing;
                private transient Stopwatch shipping;
                private long readsEmitted = 0;
                Date started;

                @Override
                public void startBundle(Context c) throws Exception {
                    super.startBundle(c);
                    if (cloudStorageUrl) {
                        storageClient = GCSOptions.Methods.createStorageClient(c.getPipelineOptions().as(GCSOptions.class), auth);
                    }
                    started = new Date();
                    readsEmitted = 0;
                    querying = Stopwatch.createUnstarted();
                    loading = Stopwatch.createUnstarted();
                    choosing = Stopwatch.createUnstarted();
                    shipping = Stopwatch.createUnstarted();
                }

                @Override
                public void finishBundle(Context c) throws Exception {
                    super.finishBundle(c);
                    final Instant instant = new Date().toInstant();
                    long ms = Duration.between(started.toInstant(), instant).toMillis();
                    String q = querying.elapsed(TimeUnit.MILLISECONDS) + " ms";
                    String l = loading.elapsed(TimeUnit.MILLISECONDS) + " ms";
                    String ch = choosing.elapsed(TimeUnit.MILLISECONDS) + " ms";
                    String s = shipping.elapsed(TimeUnit.MILLISECONDS) + " ms";
                    bunny.stepEnd("Bundle sent " + readsEmitted + " reads in " + ms + " ms. Breakdown: querying " + q + " loading " + l + " choosing " + ch + " shipping " + s);
                }

                @Override
                public void processElement(ProcessContext c) throws Exception {
                    int rejected = 0;
                    int accepted = 0;
                    Stopwatch processing = Stopwatch.createStarted();

                    Stopwatch kving = Stopwatch.createUnstarted();
                    Stopwatch outputting = Stopwatch.createUnstarted();

                    Stopwatch readingInput = Stopwatch.createStarted();
                    SimpleInterval i = c.element();
                    readingInput.stop();
                    Stopwatch init = Stopwatch.createStarted();
                    ArrayList<SimpleInterval> is = new ArrayList<>();
                    is.add(i);

                    Stopwatch cutShards = Stopwatch.createStarted();
                    List<SimpleInterval> shards = cutToShards(is, outputShardSize);
                    cutShards.stop();

                    int currentShardIndex = 0;
                    SimpleInterval currentShard = shards.get(currentShardIndex);
                    List<GATKRead> ret = new ArrayList<GATKRead>();
                    init.stop();

                    if (!cloudStorageUrl && !hadoopUrl) {
                        // local file, use setIntervalsForTraversal
                        querying.start();
                        ReadsDataSource rds = new ReadsDataSource(new File(bam));
                        rds.setIntervalsForTraversal(is);
                        querying.stop();
                        choosing.start();
                        for (GATKRead r : rds) {
                            if (!accept(r, i)) {
                                rejected++;
                                continue;
                            }
                            while (currentShard.getEnd() < r.getStart()) {
                                // move to the next shard.
                                if (!ret.isEmpty()) {
                                    accepted += ret.size();
                                    choosing.stop();
                                    shipping.start();
                                    kving.start();
                                    KV<String, Iterable<GATKRead>> out = KV.of(currentShard.toString(), ret);
                                    kving.stop();
                                    outputting.start();
                                    c.output(out);
                                    outputting.stop();
                                    //bunny.stepEnd("output subshard " + currentShardIndex + ", " + ret.size() + " reads");
                                    shipping.stop();
                                    choosing.start();
                                    ret.clear();
                                }
                                currentShardIndex++;
                                currentShard = shards.get(currentShardIndex);
                            }
                            ret.add(r);
                        }
                        choosing.stop();
                        // what remains in ret is sent after the if statement
                    } else if (cloudStorageUrl) {
                        // GCS file, load just the bit we need.
                        querying.start();
                        final SamReader reader = BAMIO.openBAM(storageClient, bam, ValidationStringency.DEFAULT_STRINGENCY);
                        SAMRecordIterator query = reader.query(i.getContig(), i.getStart(), i.getEnd(), false);
                        // we read multiple at once to estimate how long we spend waiting for disk.
                        BatchingIterator<SAMRecord> batcher = new BatchingIterator<SAMRecord>(query, outputShardSize);
                        querying.stop();

                        loading.start();
                        while (batcher.hasNext()) {
                            ArrayList<SAMRecord> buf = batcher.next();
                            loading.stop();
                            choosing.start();
                            bunny.stepEnd("done loading a batch");
                            for (SAMRecord r : buf) {
                                if (!accept(r, i)) {
                                    rejected++;
                                    continue;
                                }

                                while (currentShard.getEnd() < r.getStart()) {
                                    // move to the next shard.
                                    if (!ret.isEmpty()) {
                                        accepted += ret.size();
                                        choosing.stop();
                                        shipping.start();
                                        kving.start();
                                        KV<String, Iterable<GATKRead>> out = KV.of(currentShard.toString(), ret);
                                        kving.stop();
                                        outputting.start();
                                        c.output(out);
                                        outputting.stop();
                                        bunny.stepEnd("output subshard " + currentShardIndex + ", " + ret.size() + " reads");
                                        shipping.stop();
                                        choosing.start();
                                        ret.clear();
                                    }
                                    currentShardIndex++;
                                    currentShard = shards.get(currentShardIndex);
                                }

                                // work around some serialization madness
                                r.setHeader(null);
                                ret.add(new SAMRecordToGATKReadAdapter(r));
                            }
                            buf.clear();
                            batcher.reuse(buf);
                            choosing.stop();
                            loading.start();
                        }
                        bunny.stepEnd("that was the last batch");
                        loading.stop();
                    } else {
                        // TODO: HDFS input
                        throw new RuntimeException("Sorry, Hadoop support not yet implemented");
                    }
                    shipping.start();
                    if (!ret.isEmpty()) {
                        kving.start();
                        accepted += ret.size();
                        KV<String, Iterable<GATKRead>> out = KV.of(currentShard.toString(), ret);
                        kving.stop();
                        outputting.start();
                        c.output(out);
                        outputting.stop();
                        //bunny.stepEnd("output " + ret.size() + " reads");
                    }
                    shipping.stop();
                    readsEmitted += accepted;
                    //bunny.stepEnd("loading " + accepted + " reads for " + i.toString()
                    //    + " in " + processing.elapsed(TimeUnit.MILLISECONDS) + " ms. "
                    //    + " readingInput: " + readingInput.elapsed(TimeUnit.MILLISECONDS) + ", "
                    //    + " init: " + init.elapsed(TimeUnit.MILLISECONDS) + ", "
                    //    + " cutShards: " + cutShards.elapsed(TimeUnit.MILLISECONDS) + ", "
                    //    + " kv.of: " + kving.elapsed(TimeUnit.MILLISECONDS) + ", "
                    //    + " output: " + outputting.elapsed(TimeUnit.MILLISECONDS) + ", "
                    //    + " rejected " + rejected);
                }

                private boolean accept(GATKRead r, SimpleInterval region) {
                    if (r.isUnmapped()) {
                        // we don't want unmapped reads, *unless* their mate is in the region
                        return (!r.mateIsUnmapped() && r.getMateStart() >= region.getStart());
                    } else {
                        // the query gives us reads that overlap the region. We want only reads that *start* in those intervals.
                        return (r.getStart() >= region.getStart());
                    }
                }

                private boolean accept(SAMRecord r, SimpleInterval region) {
                    if (r.getReadUnmappedFlag() ||
                        r.getReferenceName() == null || r.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
                        r.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                        // unmapped read. We don't want it, *unless* their mate is in the region
                        return (!r.getMateUnmappedFlag() && r.getMateAlignmentStart() >= region.getStart());
                    } else {
                        // mapped read.
                        // the query gives us reads that overlap the region. We want only reads that *start* in those intervals.
                        return (r.getStart() >= region.getStart());
                    }
                }


            }));
    }


    // chr2:1-200 -> chr2:1-100,chr2:101-200
    static public List<SimpleInterval> cutToShards(Iterable<SimpleInterval> intervals, int shardSize) {
        ArrayList<SimpleInterval> ret = new ArrayList<SimpleInterval>();
        for (SimpleInterval i : intervals) {
            int beginShard = shardIndex(i.getStart(), shardSize);
            int endShard = shardIndex(i.getEnd(), shardSize);
            if (beginShard==endShard) {
                ret.add(i);
                continue;
            }
            // more than one shard: output begin to end-of-shard, then multiple full shards, then begin-of-shard to end.
            ret.add(new SimpleInterval(i.getContig(), i.getStart(), endOfShard(beginShard, shardSize)));
            for (int shard = beginShard+1; shard<endShard; shard++) {
                ret.add(new SimpleInterval(i.getContig(), beginOfShard(shard, shardSize), endOfShard(shard, shardSize)));
            }
            ret.add(new SimpleInterval(i.getContig(), beginOfShard(endShard, shardSize), i.getEnd()));
        }
        return ret;
    }

    // number of the shard this offset is in. Shards are numbered starting at zero.
    static public int shardIndex(int oneBasedOffset, int shardSize) {
        return ((oneBasedOffset-1) / shardSize);
    }

    // first offset in this shard (1-based).
    static public int beginOfShard(int shardIndex, int shardSize) {
        return ((shardIndex) * shardSize) + 1;
    }

    // last offset in this shard (1-based).
    static public int endOfShard(int shardIndex, int shardSize) {
        return beginOfShard(shardIndex+1, shardSize)-1;
    }


    /**
     * Iterator of T -> iterator of arrays of T (by reading multiple at once).
     * The last array may be smaller, as we run out of elements.
     */
    private static final class BatchingIterator<T> implements Iterator<ArrayList<T>> {
        final Iterator<T> iter;
        final int count;
        ArrayList<T> buf;

        // will give you "count" elements at a time, from "iter".
        public BatchingIterator(final Iterator<T> iter, int count)  {
            this(iter, count, null);
        }

        // if you pass it an array, it'll be used to put new content. Don't access it until
        // you get it back via next().
        public BatchingIterator(final Iterator<T> iter, int count, final ArrayList<T> buf)  {
            this.iter = iter;
            this.count = count;
            this.buf = buf;
            if (buf==null) this.buf = new ArrayList<T>(count);
            else buf.clear();
        }

        public boolean hasNext() {
            return iter.hasNext() || (null!=buf && !buf.isEmpty());
        }

        // the returned array is yours, it won't be modified from this point on.
        public ArrayList<T> next() {
            if (null==buf) buf = new ArrayList<T>();
            while (iter.hasNext() && buf.size()<count) {
                buf.add(iter.next());
            }
            ArrayList<T> ret = this.buf;
            this.buf = null;
            return ret;
        }
        // call this after next returns and before calling it again
        // if you want to provide a fresh buffer to use.
        // Don't access it until you get it back via next();
        public void reuse(ArrayList<T> buf) {
            // reject the offer if we already have data.
            if (this.buf!=null && !this.buf.isEmpty()) return;
            buf.clear();
            this.buf = buf;
        }
    }
}
