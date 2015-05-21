package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.CoderRegistry;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.coders.GenericJsonCoder;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.dataflow.utils.GenomicsDatasetOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.BQSR_Dataflow;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;

import java.io.*;
import java.nio.channels.Channels;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.List;


/**
 * First pass of the base quality score recalibration -- Generates recalibration table based on various covariates
 * (such as read group, reported quality score, machine cycle, and nucleotide context).
 * <p>
 * Dataflow version of BaseRecalibrator.
 */
@CommandLineProgramProperties(
        usage = "First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        usageShort = "Generates recalibration table",
        programGroup = ReadProgramGroup.class
)
public class DF_BaseRecalibrator extends CommandLineProgram {

    public static final String APP_NAME = "Hellbender";

    // ------------------------------------------
    // Command-line options

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private transient final BaseRecalibrationArgumentCollection BRAC = new BaseRecalibrationArgumentCollection();

    /**
     * Path to the input (BAM file). Eventually can be local or on GCS.
     */
    private String beforePath = null;


    /**
     * Path to the reference (SAM file). Eventually can be filesystem or gs://
     */
    private String referencePath = null;

    /**
     * Path to save the output to. Must be local for now.
     */
    @Argument(doc = "Path to save the output to. Must be filesystem.", optional = false)
    protected String outputTablesPath = "df-recaltables.sj";

    // ------------------------------------------
    // Dataflow-specific options

    @Argument(doc = "path to client-secrets.json", optional = false)
    protected String secretsFile;

    // not actually optional if you ask for remote exec
    @Argument(doc = "Google Cloud project name", optional = true)
    protected String project;

    // not actually optional if you ask for remote exec
    @Argument(doc = "Google Cloud Storage path for intermediate files", optional = true)
    protected String stagingLocation;

    @Argument(doc = "Number of Dataflow workers", optional = true)
    protected int numWorkers = 0;

    // ------------------------------------------

    // the inputs to BQSR
    private PCollection<Read> reads;
    private PCollection<SimpleInterval> skipIntervals;
    private SAMFileHeader header;

    /**
     * Do the work after command line has been parsed. RuntimeException may be
     * thrown by this method, and are reported appropriately.
     *
     * @return the return value or null is there is none.
     */
    @Override
    protected Object doWork() {

        try {


            // the parsing system helpfully provides a file, unaware that we might not mean a file on this machine.
            referencePath = BRAC.referenceArguments.getReferenceFile().getPath();
            // TODO: support more than one input.
            beforePath = BRAC.readArguments.getReadFiles().get(0).getPath();

            GenomicsOptions popts = makeDirectPipelineOptions(secretsFile);

            Pipeline pipeline = Pipeline.create(popts);
            CoderRegistry coderRegistry = pipeline.getCoderRegistry();
            coderRegistry.registerCoder(Read.class, GenericJsonCoder.of(Read.class));
            coderRegistry.registerCoder(SimpleInterval.class, SerializableCoder.of(SimpleInterval.class));

            // 1. prepare inputs/outputs, check arguments
            ingestLocalInputs(pipeline);
            BQSR_Dataflow.ensureReferenceIsReadable(popts, referencePath);
            BaseRecalibratorUprooted baseRecalibratorUprooted = BaseRecalibratorUprooted.fromArgs(header, BRAC);
            baseRecalibratorUprooted.checkClientArguments();

            // 2. set up computation
            PCollection<RecalibrationTables> aggregated =
                    BQSR_Dataflow.GetRecalibrationTables(header, reads, referencePath, BRAC, skipIntervals);

            if (BucketUtils.isCloudStorageUrl(outputTablesPath)) {
                saveSingleResultToGCS(aggregated, outputTablesPath);
            } else {
                // saving to a local path; this only makes sense if we're running locally.
                saveRecalibrationTables(aggregated, outputTablesPath);
            }

            // 3. Compute
            System.out.println("Computing!");
            try {
                pipeline.run();
            } catch (RuntimeException rx) {
                // Pipeline.run() helpfully wraps the underlying exception into a RuntimeException so that the stack trace shows pipeline.run at the top.
                // However, the environment here expects to see the correct exception type. So we instead throw the underlying exception, storing
                // the pipeline.run stack as "suppressed".
                Throwable t = rx.getCause();
                if (t instanceof Exception) {
                    Exception x = (Exception) t;
                    x.addSuppressed(rx);
                    throw x;
                }
                throw rx;
            }

            // 5. Get the table back and output it in text form to RAC.RECAL_TABLE.
            final InputStream tblIn;
            if (BucketUtils.isCloudStorageUrl(outputTablesPath)) {
                tblIn = Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(popts).open(GcsPath.fromUri(outputTablesPath)));
            } else {
                tblIn = new FileInputStream(outputTablesPath);
            }
            ObjectInputStream oin = new ObjectInputStream(tblIn);
            Object o = oin.readObject();
            RecalibrationTables rt = (RecalibrationTables) o;
            baseRecalibratorUprooted.onTraversalStart(null);
            baseRecalibratorUprooted.saveReport(rt, baseRecalibratorUprooted.getRequestedCovariates());


        } catch (RuntimeException rx) {
            throw rx;
        } catch (Exception x) {
            // can't change the signature.
            throw new RuntimeException(x);
        }

        return null;
    }

    // local files on the client's disk -> PCollections.
    private void ingestLocalInputs(final Pipeline pipeline) throws IOException {
        if (BucketUtils.isCloudStorageUrl(beforePath)) {
            // set up ingestion on the cloud
            // but read the header locally
            GcsPath path = GcsPath.fromUri(beforePath);
            InputStream inputstream = Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(pipeline.getOptions())
                    .open(path));
            SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(inputstream));
            header = reader.getFileHeader();
            // this won't work!
            // TODO(jpmartin): figure this out.
            reads = new ReadsSource(beforePath, pipeline).getReadPCollection(null);
        } else {
            // ingestion from local file
            SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(beforePath));
            header = reader.getFileHeader();
            List<Read> readLst = new ArrayList<>();
            ReadFilter readFilter = BaseRecalibratorUprooted.readFilter();
            for (SAMRecord sr : reader) {
                if (!readFilter.test(sr)) continue;
                try {
                    Read e = ReadConverter.makeRead(sr);
                    readLst.add(e);
                } catch (SAMException x) {
                    System.out.println("Skipping read " + sr.getReadName() + " because we can't convert it.");
                } catch (NullPointerException y) {
                    System.out.println("Skipping read " + sr.getReadName() + " because we can't convert it. (null?)");
                }
            }
            reads = pipeline.apply(Create.of(readLst).withName("input ingest"));
        }

        List<SimpleInterval> knownSitesLst = new ArrayList<>();
        for (FeatureInput<Feature> vcfSource : BRAC.RAC.knownSites) {
            FeatureDataSource<VariantContext> source = new FeatureDataSource<>(vcfSource.getFeatureFile(), new VCFCodec(), "KnownIntervals");
            for (VariantContext foo : source) {
                int start = foo.getStart();
                int end = foo.getEnd();
                String contig = foo.getContig();
                knownSitesLst.add(new SimpleInterval(contig, start, end));
            }
        }
        skipIntervals = pipeline.apply(Create.of(knownSitesLst).withName("known intervals ingest"))
                .setCoder(SerializableCoder.of(SimpleInterval.class)); // Dataflow boilerplate
    }

    private static void saveRecalibrationTables(final PCollection<RecalibrationTables> recalibrationTables, String fname) {
        recalibrationTables.apply(ParDo
                .named("saveRecalibrationTables")
                .of(new DoFn<RecalibrationTables, Void>() {
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        RecalibrationTables r = c.element();
                        System.out.println("Saving RecalibrationTables to " + fname);
                        ObjectOutputStream os = new ObjectOutputStream(new FileOutputStream(fname));
                        os.writeObject(r);
                        os.close();
                    }
                }));
    }

    private static <T> void saveSingleResult(final PCollection<T> collection, String fname) {
        collection.apply(ParDo
                .named("save to " + fname)
                .of(new DoFn<T, Void>() {
                    @Override
                    public void processElement(ProcessContext c) throws IOException {
                        T obj = c.element();
                        try (ObjectOutputStream os = new ObjectOutputStream(new FileOutputStream(fname))) {
                            os.writeObject(obj);
                        }
                    }
                }));
    }


    // a "gcsPath" starts with "gs://"
    private static <T> void saveSingleResultToGCS(final PCollection<T> collection, String gcsDestPath) {
        collection.apply(ParDo.named("save to " + gcsDestPath)
                .of(new DoFn<T, Void>() {
                    @Override
                    public void processElement(ProcessContext c) throws IOException, GeneralSecurityException {
                        GcsPath dest = GcsPath.fromUri(gcsDestPath);
                        GcsUtil gcsUtil = new GcsUtil.GcsUtilFactory().create(c.getPipelineOptions());
                        try (ObjectOutputStream out = new ObjectOutputStream(Channels.newOutputStream(gcsUtil.create(dest, "application/octet-stream")))) {
                            out.writeObject(c.element());
                        }
                    }
                }));
    }

    // Next you probably want to call: Pipeline pipeline = Pipeline.create(popts);
    private static GenomicsOptions makeDirectPipelineOptions(String secretsFilePath) throws Exception {
        PipelineOptionsFactory.register(GenomicsOptions.class);
        GenomicsOptions popts = PipelineOptionsFactory.create().as(GenomicsDatasetOptions.class);
        popts.setGenomicsSecretsFile(secretsFilePath);
        popts.setAppName(APP_NAME);
        popts.setRunner(DirectPipelineRunner.class);
        return popts;
    }


}
