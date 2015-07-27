package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.api.services.genomics.model.Read;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.dev.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.transforms.GoogleGenomicsReadToGATKRead;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.seqdoop.hadoop_bam.util.SAMHeaderReader;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static java.util.Arrays.stream;

/**
 * Class to load reads into a PCollection from a cloud storage bucket, a Hadoop filesystem, or a local bam file.
 */
public final class ReadsDataflowSource implements Serializable {
    private final String bam;
    private final boolean cloudStorageUrl;
    private final boolean hadoopUrl;
    private transient GCSOptions options;
    private transient Pipeline pipeline;
    private GenomicsFactory.OfflineAuth auth;
    private final static Logger logger = LogManager.getLogger(ReadsDataflowSource.class);
    private static final long serialVersionUID = 1L;

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
        } catch (IOException e) {
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
        return getReadPCollection(intervals, ValidationStringency.SILENT);
    }

    /**
     * Create a {@link PCollection<GATKRead>} containing all the reads overlapping the given intervals.
     * Reads that are unmapped are ignored.
     * @param intervals a list of SimpleIntervals.  These must be non-overlapping intervals or the results are undefined.
     * @param stringency how to react to malformed reads.
     * @return a PCollection containing all the reads that overlap the given intervals.
     */
    public PCollection<GATKRead> getReadPCollection(List<SimpleInterval> intervals, ValidationStringency stringency) {
        PCollection<GATKRead> preads;
        if(cloudStorageUrl){
            Iterable<Contig> contigs = intervals.stream()
                    .map(i -> new Contig(i.getContig(), i.getStart(), i.getEnd()))
                    .collect(Collectors.toList());

            PCollection<Read> rawReads = ReadBAMTransform.getReadsFromBAMFilesSharded(pipeline, auth, contigs, stringency, ImmutableList.of(bam));
            preads = rawReads.apply(new GoogleGenomicsReadToGATKRead());
        } else if (hadoopUrl) {
            preads = DataflowUtils.getReadsFromHadoopBam(pipeline, intervals, stringency, bam);
        } else {
            preads = DataflowUtils.getReadsFromLocalBams(pipeline, intervals, stringency, ImmutableList.of(new File(bam)));
        }
        return preads;
    }

    public static PCollectionView<SAMFileHeader> getHeaderView(Pipeline p, SAMFileHeader readsHeader) {
        return p.apply(Create.of(readsHeader).setName("reads header")).setCoder(SerializableCoder.of(SAMFileHeader.class)).apply(View.<SAMFileHeader>asSingleton());
    }


    /**
     * Create a {@link PCollection<KV<String,Iterable<GATKRead>>>} containing all the reads starting in the given intervals,
     * grouped by the shard they start in.
     * The shards are strings as per SimpleInterval's ToString syntax.
     * Reads that are malformed or unmapped are ignored.
     * @param intervals a list of SimpleIntervals.
     */
    public PCollection<KV<String,Iterable<GATKRead>>> getGroupedReadPCollection(List<SimpleInterval> intervals, int shardSize, Pipeline pipeline) {
        List<SimpleInterval> shards = cutToShards(intervals, shardSize);
        logger.info("Input range sharded into "+shards.size()+" intervals.");
        List<KV<String,Void>> intervalKeys;
        intervalKeys = shards.stream().map(s-> KV.<String,Void>of(s.toString(), null)).collect(Collectors.<KV<String,Void>>toList());
        return pipeline
            // ship shards to workers
            .apply(Create.<KV<String, Void>>of(intervalKeys).setName("shards"))
            // group by shard
            .apply(GroupByKey.create())
            // for each shard, load the corresponding reads
            .apply(ParDo.named("getGroupedReadPCollection").of(new DoFnWLog<KV<String, Iterable<Void>>, KV<String, Iterable<GATKRead>>>("getGroupedReadPCollection") {
                private static final long serialVersionUID = 1L;

                @Override
                public void processElement(ProcessContext c) throws Exception {
                    int rejected = 0;
                    String intStr = c.element().getKey();
                    SimpleInterval i = new SimpleInterval(intStr);
                    List<GATKRead> ret = new ArrayList<GATKRead>();
                    Storage.Objects storageClient = GCSOptions.Methods.createStorageClient(c.getPipelineOptions().as(GCSOptions.class), auth);
                    final SamReader reader = BAMIO.openBAM(storageClient, bam, ValidationStringency.DEFAULT_STRINGENCY);
                    SAMRecordIterator query = reader.query(i.getContig(), i.getStart(), i.getEnd(), false);
                    while (query.hasNext()) {
                        SAMRecord r = query.next();
                        // the query gives us reads that overlap the region. We want only reads that *start* in those intervals.
                        if (r.getStart() < i.getStart()) {
                            rejected++;
                            continue;
                        }
                        // we also don't want unmapped reads
                        if (r.getReadUnmappedFlag() ||
                            r.getReferenceName() == null || r.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
                            r.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                            rejected++;
                            continue;
                        }
                        Read genR = ReadConverter.makeRead(r);
                        GATKRead gatkR = new GoogleGenomicsReadToGATKReadAdapter(genR);
                        if (gatkR.isUnmapped()) {
                            rejected++;
                            continue;
                        }
                        ret.add(gatkR);
                    }
                    c.output(KV.of(intStr, ret));
                    bunny.stepEnd("loading " + ret.size() + " reads for " + intStr + ", rejected " + rejected);
                }
            }));
    }

    // chr2:1-200 -> chr2:1-100,chr2:101-200
    static public List<SimpleInterval> cutToShards(List<SimpleInterval> intervals, int shardSize) {
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
}