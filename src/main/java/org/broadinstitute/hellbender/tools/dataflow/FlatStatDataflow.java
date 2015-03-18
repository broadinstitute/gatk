package org.broadinstitute.hellbender.tools.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.SearchReadsRequest;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.options.*;
import com.google.cloud.dataflow.sdk.transforms.Count;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.ReadReader;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.readers.bam.Reader;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsDatasetOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.cloud.genomics.utils.Paginator;
import com.google.common.base.Function;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.IntervalArgumentCollection;
import org.codehaus.jackson.annotate.JsonIgnore;

import java.io.IOException;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static com.google.common.collect.Lists.newArrayList;

public class FlatStatDataflow {

    private static final Logger LOG = Logger.getLogger(FlatStatDataflow.class.getName());
    private static CountReadsOptions options;
    private static Pipeline p;
    private static GenomicsFactory.OfflineAuth auth;

    @Argument(doc="gs:// path to a bam file")
    public String bamPath;

    @ArgumentCollection
    public IntervalArgumentCollection intervalCollection = new IntervalArgumentCollection();

    public static interface CountReadsOptions extends GenomicsDatasetOptions, GCSOptions {
        @Description("The ID of the Google Genomics ReadGroupSet this pipeline is working with. "
                + "Default (empty) indicates all ReadGroupSets.")
        @Default.String("")
        String getReadGroupSetId();

        void setReadGroupSetId(String readGroupSetId);

        @Description("The path to the BAM file to get reads data from.")
        @Default.String("")
        String getBAMFilePath();

        void setBAMFilePath(String filePath);

        @Description("Whether to shard BAM reading")
        @Default.Boolean(true)
        boolean getShardBAMReading();

        void setShardBAMReading(boolean newValue);

        class ContigsFactory implements DefaultValueFactory<Iterable<Contig>> {
            @Override
            public Iterable<Contig> create(PipelineOptions options) {
                return Iterables.transform(Splitter.on(",").split(options.as(CountReadsOptions.class).getReferences()),
                        new Function<String, Contig>() {
                            @Override
                            public Contig apply(String contigString) {
                                ArrayList<String> contigInfo = newArrayList(Splitter.on(":").split(contigString));
                                return new Contig(contigInfo.get(0),
                                        contigInfo.size() > 1 ?
                                                Long.valueOf(contigInfo.get(1)) : 0,
                                        contigInfo.size() > 2 ?
                                                Long.valueOf(contigInfo.get(2)) : -1);
                            }
                        });
            }
        }

        @Default.InstanceFactory(ContigsFactory.class)
        @JsonIgnore
        Iterable<Contig> getContigs();

        void setContigs(Iterable<Contig> contigs);
    }

    public static void main(String[] args) throws GeneralSecurityException, IOException {
        options = PipelineOptionsFactory.fromArgs(args).withValidation().as(CountReadsOptions.class);
        GenomicsOptions.Methods.validateOptions(options);
        auth = GCSOptions.Methods.createGCSAuth(options);
        p = Pipeline.create(options);
        DataflowWorkarounds.registerGenomicsCoders(p);

        PCollection<Read> reads = getReads();
        PCollection<Long> readCount = reads.apply(Count.<Read>globally());
        PCollection<String> readCountText = readCount.apply(ParDo.of(new DoFn<Long, String>() {
            @Override
            public void processElement(DoFn<Long, String>.ProcessContext c) throws Exception {
                c.output(String.valueOf(c.element()));
            }
        }).named("toString"));
        readCountText.apply(TextIO.Write.to(options.getOutput()).named("WriteOutput"));
        p.run();
    }

    private static PCollection<Read> getReads() throws IOException {
        if (!options.getBAMFilePath().isEmpty()) {
            return getReadsFromBAMFile();
        }
        if (!options.getReadGroupSetId().isEmpty()) {
            return getReadsFromAPI();
        }
        throw new IOException("Either BAM file or ReadGroupSet must be specified");
    }

    private static PCollection<Read> getReadsFromAPI() {
        List<SearchReadsRequest> requests = getReadRequests(options);
        PCollection<SearchReadsRequest> readRequests =
                DataflowWorkarounds.getPCollection(requests, p, options.getNumWorkers());
        PCollection<Read> reads =
                readRequests.apply(
                        ParDo.of(
                                new ReadReader(auth, Paginator.ShardBoundary.OVERLAPS))
                                .named(ReadReader.class.getSimpleName()));
        return reads;
    }

    private static List<SearchReadsRequest> getReadRequests(CountReadsOptions options) {
        List<SearchReadsRequest> requests = Lists.newArrayList();
        requests.add(new SearchReadsRequest()
                .setReadGroupSetIds(
                        Collections.singletonList(options.getReadGroupSetId()))
                .setReferenceName(options.getReferences())
                .setPageSize(2048));

        return requests;
    }

    private static PCollection<Read> getReadsFromBAMFile() throws IOException {
        LOG.info("getReadsFromBAMFile");

        final Iterable<Contig> contigs = options.getContigs();

        if (options.getShardBAMReading()) {
            return ReadBAMTransform.getReadsFromBAMFilesSharded(p,
                    auth,
                    contigs,
                    Collections.singletonList(options.getBAMFilePath()));
        } else {  // For testing and comparing sharded vs. not sharded only
            return p.apply(
                    Create.of(
                            Reader.readSequentiallyForTesting(
                                    GCSOptions.Methods.createStorageClient(options, auth),
                                    options.getBAMFilePath(),
                                    contigs.iterator().next())));
        }
    }
}
