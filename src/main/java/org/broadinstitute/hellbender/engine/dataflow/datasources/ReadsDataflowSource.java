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
import org.broadinstitute.hellbender.engine.dataflow.transforms.GoogleGenomicsReadToGATKRead;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.seqdoop.hadoop_bam.util.SAMHeaderReader;

import java.io.File;
import java.io.IOException;
import java.security.GeneralSecurityException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class to load reads into a PCollection from a cloud storage bucket, a Hadoop filesystem, or a local bam file.
 */
public final class ReadsDataflowSource {
    private final String bam;
    private final boolean cloudStorageUrl;
    private final boolean hadoopUrl;
    private GCSOptions options;
    private final Pipeline pipeline;
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
        return ShardedReadsDatasource.get(bam, intervals, readShardSize, outputShardSize, pipeline, auth);
    }

}
