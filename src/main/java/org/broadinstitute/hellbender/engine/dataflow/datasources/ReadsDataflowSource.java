package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.api.services.genomics.model.Read;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.readers.bam.ReaderOptions;
import com.google.cloud.genomics.dataflow.readers.bam.ShardingPolicy;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.base.Stopwatch;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.dataflow.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.transforms.GoogleGenomicsReadToGATKRead;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.seqdoop.hadoop_bam.util.SAMHeaderReader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.security.GeneralSecurityException;
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
 * Class to load reads into a PCollection from a cloud storage bucket, a Hadoop filesystem, or a local bam file.
 */
public final class ReadsDataflowSource implements Serializable {
    private static final long serialVersionUID = 1l;
    private final static Logger logger = LogManager.getLogger(ReadsDataflowSource.class);
    private final String bam;
    private final boolean cloudStorageUrl;
    private final boolean hadoopUrl;
    private transient GCSOptions options;
    private transient final Pipeline pipeline;
    private GenomicsFactory.OfflineAuth auth;

    /**
     * @param bam A file path or a google bucket identifier to a bam file to read
     * @param p the pipeline object for the job. This is needed to read a bam from a bucket.
     *          The options inside of the pipeline MUST BE GCSOptions (to get the secret file).
     */
    public ReadsDataflowSource(String bam, Pipeline p){
        this.bam = Utils.nonNull(bam);
        this.pipeline = p;

        cloudStorageUrl = BucketUtils.isCloudStorageUrl(bam);
        hadoopUrl = BucketUtils.isHadoopUrl(bam);
        if(cloudStorageUrl) {
            // The options used to create the pipeline must be GCSOptions to get the secret file.
            try {
                options = p.getOptions().as(GCSOptions.class);
            } catch (ClassCastException e) {
                throw new GATKException("The pipeline options was not GCSOptions.", e);
            }
            GenomicsOptions.Methods.validateOptions(options);
            auth = getAuth(options);
        }
    }

    private static GenomicsFactory.OfflineAuth getAuth(GCSOptions options){
        try {
            return GCSOptions.Methods.createGCSAuth(options);
        } catch (IOException | GeneralSecurityException e) {
            throw new GATKException("Couldn't create a dataflow auth object.", e);
        }
    }

    /**
     * Get a SAMFileHeader to use with Reads produced by this ReadsSource
     */
    public SAMFileHeader getHeader() {
        if(cloudStorageUrl) {
            try {
                Storage.Objects storageClient = GCSOptions.Methods.createStorageClient(options, auth);
                final SamReader reader = BAMIO.openBAM(storageClient, bam, ValidationStringency.DEFAULT_STRINGENCY);
                return reader.getFileHeader();
            } catch (IOException e) {
                throw new GATKException("Failed to read bams header from " + bam + ".", e);
            }
        } else if (hadoopUrl) {
            try {
                return SAMHeaderReader.readSAMHeaderFrom(new Path(bam), new Configuration());
            } catch (IOException e) {
                throw new GATKException("Failed to read bams header from " + bam + ".", e);
            }
        } else {
            return SamReaderFactory.makeDefault().getFileHeader(new File(bam));
        }
    }

    /**
     * Create a {@link PCollection<GATKRead>} containing all the reads overlapping the given intervals.
     * Reads that are malformed or unmapped are ignored.
     * @param intervals a list of SimpleIntervals.  These must be non-overlapping intervals or the results are undefined.
     * @return a PCollection containing all the reads that overlap the given intervals.
     */
    public PCollection<GATKRead> getReadPCollection(List<SimpleInterval> intervals) {
        return getReadPCollection(intervals, ValidationStringency.SILENT, false);
    }

    /**
     * Create a {@link PCollection<GATKRead>} containing all the reads overlapping the given intervals.
     * Reads that are unmapped are ignored.
     * @param intervals a list of SimpleIntervals.  These must be non-overlapping intervals or the results are undefined.
     * @param stringency how to react to malformed reads.
     * @param includeUnmappedReads to include unmapped reads.
     * @return a PCollection containing all the reads that overlap the given intervals.
     */
    public PCollection<GATKRead> getReadPCollection(List<SimpleInterval> intervals, ValidationStringency stringency, boolean includeUnmappedReads) {
        PCollection<GATKRead> preads;
        if(cloudStorageUrl){
            Iterable<Contig> contigs = intervals.stream()
                    .map(i -> new Contig(i.getContig(), i.getStart(), i.getEnd()))
                    .collect(Collectors.toList());
            try {
                PCollection<Read> rawReads = ReadBAMTransform.getReadsFromBAMFilesSharded(pipeline, auth,contigs, new ReaderOptions(stringency, includeUnmappedReads), bam, ShardingPolicy.LOCI_SIZE_POLICY);
                preads = rawReads.apply(new GoogleGenomicsReadToGATKRead());
            } catch (IOException ex) {
                throw new UserException.CouldNotReadInputFile("Unable to read "+bam, ex);
            }
        } else if (hadoopUrl) {
            preads = DataflowUtils.getReadsFromHadoopBam(pipeline, intervals, stringency, bam);
        } else {
            preads = DataflowUtils.getReadsFromLocalBams(pipeline, intervals, stringency, ImmutableList.of(new File(bam)));
        }
        return preads;
    }

    public static PCollectionView<SAMFileHeader> getHeaderView(Pipeline p, SAMFileHeader readsHeader) {
        return p.apply(Create.of(readsHeader).withCoder(SerializableCoder.of(SAMFileHeader.class))).setName("reads header").apply(View.<SAMFileHeader>asSingleton());
    }



    /**
     * Create a {@link PCollection< KV <String,Iterable<GATKRead>>>} containing all the reads starting in the given intervals and their mates,
     * grouped by the shard they start in.
     * The shards are strings as per SimpleInterval's ToString syntax.
     * Reads that are malformed or unmapped are ignored.
     *
     * Note: for performance reasons, this will give you SAMRecords that have no header. You have to put it back in.
     *
     * @param intervals a list of SimpleIntervals.
     * @param readShardSize shards of this size are assigned to individual machines (streamed in from GCS).
     * @param  outputShardSize shards of this many reads are output at a time (held all together in RAM).
     */
    public PCollection<KV<String,Iterable<GATKRead>>> getGroupedReadPCollection(List<SimpleInterval> intervals, int readShardSize, int outputShardSize, Pipeline pipeline) {
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

                                //Read genR = ReadConverter.makeRead(r);
                                //GATKRead gatkR = new GoogleGenomicsReadToGATKReadAdapter(genR);
                                //ret.add(gatkR);

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


    public class BatchingIterator<T> implements Iterator<ArrayList<T>> {
        Iterator<T> iter; int count; ArrayList<T> buf;
        public BatchingIterator(Iterator<T> iter, int count)  {
            this(iter, count, null);
        }
        public BatchingIterator(Iterator<T> iter, int count, ArrayList<T> buf)  {
            this.iter = iter;
            this.count = count;
            this.buf = buf;
            if (buf==null) this.buf = new ArrayList<T>(count);
        }
        public boolean hasNext() {
            return iter.hasNext() || (null!=buf && !buf.isEmpty());
        }
        public ArrayList<T> next() {
            if (null==buf) buf = new ArrayList<T>();
            while (iter.hasNext() && buf.size()<count) {
                buf.add(iter.next());
            }
            ArrayList<T> ret = this.buf;
            this.buf = null;
            return ret;
        }
        public void reuse(ArrayList<T> buf) {
            // reject the offer if we already have data.
            if (this.buf!=null && !this.buf.isEmpty()) return;
            buf.clear();
            this.buf = buf;
        }
    }
}
