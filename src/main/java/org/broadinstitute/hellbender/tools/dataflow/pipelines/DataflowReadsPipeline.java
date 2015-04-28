package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.transforms.Filter;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowSAMFn;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.dataflow.SAMSerializableFunction;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Handles lifting reads from a bams into a PCollection and provides hooks to apply {@link ReadFilter}, {@link ReadConverter}, and
 * a {@link PTransformSAM}.
 *
 *Subclasses must override {@link #getTool()} and optionally override {@link #getReadFilters()} and {@link #getReadTransformers()}
 */
public abstract class DataflowReadsPipeline extends DataflowCommandLineProgram {

    @Argument(doc="a prefix for the dataflow output files", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "path to the client secret file for google cloud authentication, necessary if accessing data from buckets",
            shortName = "secret", fullName = "client_secret", optional=true)
    protected File clientSecret = new File("client_secret.json");

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new IntervalArgumentCollection();

    /**
     * Returns a transform which performs the work of the given tool.
     * Subclasses must override this.
     * @return the Transform that should be performed in this pipeline.
     */
    abstract protected PTransformSAM<?> getTool();

    /**
     * @return {@link ReadFilter}s to apply before running the tool
     */
    protected ImmutableList<ReadFilter> getReadFilters(){
        return ImmutableList.of();
    }

    /**
     * @return {@link ReadTransformer}s to apply before running the tool
     */
    protected ImmutableList<ReadTransformer> getReadTransformers(){
        return ImmutableList.of();
    }

    private static SerializableFunction<Read,Boolean> wrapFilter(ReadFilter filter, String headerString){
        return new SAMSerializableFunction<>(headerString, filter);
    }

    @Override
    final protected void setupPipeline(Pipeline pipeline) {
        ReadsSource readsSource = new ReadsSource(bam, clientSecret);
        String headerString = readsSource.getHeaderString();
        List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(ReadUtils.samHeaderFromString(headerString).getSequenceDictionary());

        PCollection<Read> preads = readsSource.getReadPCollection(pipeline, intervals);

        PCollection<?> presult = applyTransformsToPipeline(headerString, preads);

        PCollection<String> pstrings = presult.apply(DataflowUtils.convertToString());
        pstrings.apply(TextIO.Write.to(outputFile));
    }

    @VisibleForTesting
    protected PCollection<?> applyTransformsToPipeline(String headerString, PCollection<Read> preads) {
        for (ReadFilter filter : getReadFilters()) {
            preads = preads.apply(Filter.by(wrapFilter(filter, headerString)));
        }

        for (ReadTransformer transformer : getReadTransformers()){
            preads = preads.apply(wrapTransformer(transformer, headerString));
        }

        PTransformSAM<?> f = getTool();
        f.setHeaderString(headerString);

        return preads.apply(f);
    }

    private static PTransform<? super PCollection<Read>,PCollection<Read>> wrapTransformer(ReadTransformer transformer, String headerString){
        return ParDo.of(new DataFlowSAMFn<Read>(headerString){

                    @Override
                    protected void apply(SAMRecord read) {
                        output(ReadConverter.makeRead(transformer.apply(read)));
                    }
                }
        );
    }

    private static class ReadsSource{
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

}
