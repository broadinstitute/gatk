package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.seqdoop.hadoop_bam.util.SAMHeaderReader;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class to load reads into a PCollection from a cloud storage bucket, a Hadoop filesystem, or a local bam file.
 */
public final class ReadsSource {
    private final String bam;
    private final boolean cloudStorageUrl;
    private final boolean hadoopUrl;
    private GCSOptions options;
    private Pipeline pipeline;
    private GenomicsFactory.OfflineAuth auth;

    /**
     * @param bam A file path or a google bucket identifier to a bam file to read
     * @param p the pipeline object for the job. This is needed to read a bam from a bucket.
     *          The options inside of the pipeline MUST BE GCSOptions (to get the secret file).
     */
    public ReadsSource(String bam, Pipeline p){
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
     * Create a {@link PCollection<Read>} containing all the reads overlapping the given intervals.
     * Reads that are malformed or unmapped are ignored.
     * @param intervals a list of SimpleIntervals.  These must be non-overlapping intervals or the results are undefined.
     * @return a PCollection containing all the reads that overlap the given intervals.
     */
    public PCollection<Read> getReadPCollection(List<SimpleInterval> intervals) {
        return getReadPCollection(intervals, ValidationStringency.SILENT);
    }

    /**
     * Create a {@link PCollection<Read>} containing all the reads overlapping the given intervals.
     * Reads that are unmapped are ignored.
     * @param intervals a list of SimpleIntervals.  These must be non-overlapping intervals or the results are undefined.
     * @param stringency how to react to malformed reads.
     * @return a PCollection containing all the reads that overlap the given intervals.
     */
    public PCollection<Read> getReadPCollection(List<SimpleInterval> intervals, ValidationStringency stringency) {
        PCollection<Read> preads;
        if(cloudStorageUrl){
            Iterable<Contig> contigs = intervals.stream()
                    .map(i -> new Contig(i.getContig(), i.getStart(), i.getEnd()))
                    .collect(Collectors.toList());

            preads = ReadBAMTransform.getReadsFromBAMFilesSharded(pipeline, auth, contigs, stringency, ImmutableList.of(bam));
        } else if (hadoopUrl) {
            preads = DataflowUtils.getReadsFromHadoopBam(pipeline, intervals, stringency, bam);
        } else {
            preads = DataflowUtils.getReadsFromLocalBams(pipeline, intervals, stringency, ImmutableList.of(new File(bam)));
        }
        return preads;
    }
}