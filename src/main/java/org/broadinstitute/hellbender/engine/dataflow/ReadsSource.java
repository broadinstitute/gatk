package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class to load reads into a PCollection from either a cloud storage bucket or a local bam file.
 */
public class ReadsSource {
    private final String bam;
    private final boolean cloudStorageUrl;
    private GCSOptions options;
    private GenomicsFactory.OfflineAuth auth;

    public ReadsSource(String bam, File clientSecret){
        this.bam = bam;

        cloudStorageUrl = BucketUtils.isCloudStorageUrl(bam);
        if(cloudStorageUrl) {
            if (!clientSecret.exists()) {
                throw new UserException("You must specify a valid client secret file if using bams from a google bucket");
            }
            //HACK this is gross but it seemed like the easiest way to deal with the auth stuff
            options = PipelineOptionsFactory.fromArgs(new String[]{"--genomicsSecretsFile=" + clientSecret.getAbsolutePath()}).as(GCSOptions.class);
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

    public String getHeaderString() {
        if(cloudStorageUrl) {
            try {
                Storage.Objects storageClient = GCSOptions.Methods.createStorageClient(options, auth);
                final SamReader reader = BAMIO.openBAM(storageClient, bam);
                return reader.getFileHeader().getTextHeader();
            } catch (IOException e) {
                throw new GATKException("Failed to read bams header from " + bam + ".", e);
            }
        } else {
            return SamReaderFactory.makeDefault().getFileHeader(new File(bam)).getTextHeader();
        }
    }

    public PCollection<Read> getReadPCollection(Pipeline pipeline, List<SimpleInterval> intervals) {
        PCollection<Read> preads;
        if(cloudStorageUrl){
            Iterable<Contig> contigs = intervals.stream()
                    .map(i -> new Contig(i.getContig(), i.getStart(), i.getEnd()))
                    .collect(Collectors.toList());

            preads = ReadBAMTransform.getReadsFromBAMFilesSharded(pipeline, auth, contigs, ImmutableList.of(bam));
        } else {
            preads = DataflowUtils.getReadsFromLocalBams(pipeline, intervals, ImmutableList.of(new File(bam)));
        }
        return preads;
    }
}