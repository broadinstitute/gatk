package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.dev.BunnyLog;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.BaseRecalibratorDataflowUtils;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.DataflowReadFilter;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.nio.channels.Channels;
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
public class BaseRecalibratorDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    private final static Logger logger = LogManager.getLogger(BaseRecalibratorDataflow.class);
    // temporary file with the serialized recalibrationTables.
    private final static String TEMP_RECALTABLES = "temp-ds-recaltables";

    // ------------------------------------------
    // Command-line options

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private transient final BaseRecalibrationArgumentCollection BRAC = new BaseRecalibrationArgumentCollection();

    /**
     * Path to save the serialized tables to. Local or GCS.
     */
    @Argument(doc = "Path to save the serialized recalibrationTables to. If running on the cloud, either leave outputTablesPath unset or point it to a GCS location.",
            shortName = "otp", fullName = "outputTablesPath", optional = true)
    protected String outputTablesPath = null;

    // ------------------------------------------

    /**
     * Path to the reference (SAM file). Eventually can be filesystem or gs://
     */
    private String referencePath = null;

    // Whether we run on the cloud.
    private boolean remote = false;
    // Whether we want to save the textual version of the output.
    private boolean saveTextualTables;

    // the inputs to BQSR
    private PCollection<GATKRead> reads;
    private PCollection<SimpleInterval> skipIntervals;
    private SAMFileHeader header;

    // set up with the user's arguments, kept around to save the textual report.
    private BaseRecalibratorWorker baseRecalibratorWorker;

    private BunnyLog bunny = new BunnyLog(logger);

    @Override
    protected void setupPipeline(Pipeline pipeline) {
        try {
            bunny.start("BaseRecalibratorDataflow");
            saveTextualTables = (null != BRAC.RAC.RECAL_TABLE_FILE);

            referencePath = BRAC.referenceArguments.getReferenceFileName();
            reads = ingestReadsAndGrabHeader(pipeline, BRAC.readArguments.getReadFilesNames());
            skipIntervals = ingestKnownIntervals(pipeline, BRAC.RAC.knownSites);
            BaseRecalibratorDataflowUtils.ensureReferenceIsReadable(pipeline.getOptions(), referencePath);
            baseRecalibratorWorker = BaseRecalibratorWorker.fromArgs(header, BRAC);
            baseRecalibratorWorker.checkClientArguments();

            // 2. set up computation
            PCollection<RecalibrationTables> aggregated =
                    BaseRecalibratorDataflowUtils.getRecalibrationTables(header, reads, referencePath, BRAC, skipIntervals);

            // If saving textual output then we need to make sure we can get to the output
            if (saveTextualTables) {
                if (null == outputTablesPath) {
                    // we need those, so let's pick a location for them.
                    outputTablesPath = pickOutputTablesPath(isRemote(), stagingLocation);
                }
                DataflowUtils.SaveDestination dest = DataflowUtils.serializeSingleObject(aggregated, outputTablesPath);
                if (isRemote() && dest == DataflowUtils.SaveDestination.LOCAL_DISK) {
                    throw new UserException("If running on the cloud, either leave outputTablesPath unset or point it to a GCS location.");
                }
            }
            bunny.stepEnd("setup");
        } catch (UserException rx) {
            throw rx;
        } catch (GATKException rx) {
            throw rx;
        } catch (Exception x) {
            throw new GATKException("Unexpected: " + x.getMessage(), x);
        }
    }

    @Override
    protected void afterPipeline(Pipeline p) {
        bunny.stepEnd("dataflow");
        if (saveTextualTables) {
            //  Get the table back and output it in text form to RAC.RECAL_TABLE.
            // TODO: if running on the cloud and the output destination is on the cloud, then it's faster to have a worker do it directly, without the file roundtrip.
            try (ObjectInputStream oin = new ObjectInputStream(BucketUtils.openFile(outputTablesPath, p.getOptions()))) {
                Object o = oin.readObject();
                RecalibrationTables rt = (RecalibrationTables) o;
                baseRecalibratorWorker.onTraversalStart(null);
                baseRecalibratorWorker.saveReport(rt, baseRecalibratorWorker.getRequestedCovariates());
                bunny.stepEnd("repatriate_report");
            } catch (Exception e) {
                throw new GATKException("Unexpected: unable to read results file. (bug?)", e);
            }
        }
        bunny.end();
    }

    protected static String pickOutputTablesPath(boolean remote, String stagingLocation) {
        if (remote) {
            return GcsPath.fromUri(stagingLocation).resolve(TEMP_RECALTABLES + ".sj").toString();
        } else {
            File outputTables = BaseTest.createTempFile(TEMP_RECALTABLES, ".sj");
            return outputTables.getPath();
        }
    }

    /**
     * reads local disks or GCS -> header, and PCollection
     */
    private PCollection<GATKRead> ingestReadsAndGrabHeader(final Pipeline pipeline, List<String> filenames) throws IOException {
        if (filenames.size() > 1) {
            throw new UserException("Sorry, we only support a single input file for now.");
        }
        // TODO: support more than one input.
        String beforePath = filenames.get(0);

        // input reads
        if (BucketUtils.isCloudStorageUrl(beforePath)) {
            // set up ingestion on the cloud
            // but read the header locally
            GcsPath path = GcsPath.fromUri(beforePath);
            InputStream inputstream = Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(pipeline.getOptions())
                    .open(path));
            SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(inputstream));
            header = reader.getFileHeader();

            final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
            final ReadFilter readFilter = BaseRecalibratorWorker.readFilter(header);
            final List<SimpleInterval> intervals = BRAC.intervalArgumentCollection.intervalsSpecified() ? BRAC.intervalArgumentCollection.getIntervals(sequenceDictionary) :
                    IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
            return new ReadsSource(beforePath, pipeline).getReadPCollection(intervals, ValidationStringency.SILENT)
                    // keep only the ones BQSR's interested in.
                    .apply(new DataflowReadFilter(readFilter, header));
        } else {
            // ingestion from local file
            try ( ReadsDataSource readsSource = new ReadsDataSource(new File(beforePath)) ) {
                header = readsSource.getHeader();
                List<GATKRead> readLst = new ArrayList<>();
                ReadFilter readFilter = BaseRecalibratorWorker.readFilter(header);
                for ( GATKRead read : readsSource ) {
                    if ( readFilter.test(read) ) {
                        readLst.add(read);
                    }
                }
                return pipeline.apply(Create.of(readLst).setName("input ingest")).setCoder(new GATKReadCoder());
            }
        }
    }

    /** list of known intervals -> PCollection */
    @SuppressWarnings("unchecked")
    private static PCollection<SimpleInterval> ingestKnownIntervals(final Pipeline pipeline, List<FeatureInput<Feature>> knownSites) {
        // known sites
        List<SimpleInterval> knownSitesLst = new ArrayList<>();
        for (FeatureInput<Feature> vcfSource : knownSites) {
            File featureFile = vcfSource.getFeatureFile();
            FeatureDataSource<Feature> source = new FeatureDataSource<Feature>(featureFile, (FeatureCodec<Feature, ?>) FeatureManager.getCodecForFile(featureFile), "KnownIntervals");
            for (Feature f : source) {
                knownSitesLst.add(new SimpleInterval(f));
            }
        }
        return pipeline.apply(Create.of(knownSitesLst).setName("known intervals ingest"))
                .setCoder(SerializableCoder.of(SimpleInterval.class)); // Dataflow boilerplate
    }

}
