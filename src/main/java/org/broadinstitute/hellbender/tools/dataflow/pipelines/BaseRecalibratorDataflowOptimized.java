package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantsDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.transforms.composite.AddContextDataToRead;
import org.broadinstitute.hellbender.engine.dataflow.transforms.composite.AddContextDataToReadOptimized;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorFn;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorOptimizedTransform;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.logging.BunnyLog;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.File;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


/**
 * Optimized version of the base quality score recalibration -- Generates recalibration table based on various covariates
 * (such as read group, reported quality score, machine cycle, and nucleotide context).
 * <p>
 * Note: unmapped reads are ignored by this version of the tool.
 * <p>
 * This version uses the AddContextDataToRead transform.
 */
@CommandLineProgramProperties(
        summary = "First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        oneLineSummary = "Generates recalibration table",
        programGroup = ReadProgramGroup.class
)
public class BaseRecalibratorDataflowOptimized extends DataflowCommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
      * Reference window function for BQSR. For each read, returns an interval representing the span of
      * reference bases required by the BQSR algorithm for that read. Should be passed into the
      * {@link ReferenceDataflowSource} object for the {@link AddContextDataToRead} transform.
      */
    public static final SerializableFunction<GATKRead, SimpleInterval> BQSR_REFERENCE_WINDOW_FUNCTION =
            read -> BAQ.getReferenceWindowForRead(read, BAQ.DEFAULT_BANDWIDTH, BAQ.DEFAULT_INCLUDE_CLIPPED_BASES);

    public static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";

    private final static Logger logger = LogManager.getLogger(BaseRecalibratorDataflowOptimized.class);
    // temporary file with the serialized recalibrationTables.
    private final static String TEMP_RECALTABLES = "temp-ds-recaltables";

    // the granularity at which we'll want to assign work and read inputs.
    private final int bigShardSize = 1_000_000;
    // the granularity at which we want to batch reads,variants,and reference for processing.
    private final int smallShardSize = 5_000;
    // we make the assumption that all reads that start in a shard
    // will end within "margin" of the shard. This allows us to split the computation
    // across machines.
    private final int margin = 1000;


    // ------------------------------------------
    // Command-line options

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection recalArgs = new RecalibrationArgumentCollection();

    @ArgumentCollection
    private final RequiredReadInputArgumentCollection readArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    private final IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @ArgumentCollection
    private final ReferenceInputArgumentCollection referenceArguments = new RequiredReferenceInputArgumentCollection();

    @Argument(doc = "the known variants", shortName = "knownSites", fullName = "knownSites", optional = false)
    private List<String> knownVariants;

    /**
     * Path to save the temporary serialized recalibration tables to. Local or GCS.
     */
    @Argument(doc = "Path to save the serialized recalibrationTables to. If running on the cloud, either leave serializedOutputTablesPath unset or point it to a GCS location.",
            shortName = "sotp", fullName = "serializedOutputTablesPath", optional = true)
    protected String serializedOutputTablesPath = null;

    /**
     * Path to save the final recalibration tables to. Local.
     */
    @Argument(doc = "Path to save the final recalibrationTables to.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputTablesPath = null;

    // ------------------------------------------

    private SAMFileHeader readsHeader;

    private transient BunnyLog bunny = new BunnyLog(logger);

    @Override
    protected String jobName() {
        return "BaseRecalibratorDataflowOptimized-"+ (100+new Random().nextInt(900));
    }

    @Override
    protected void setupPipeline(Pipeline pipeline) {
        try {
            bunny.start("BaseRecalibratorDataflowOptimized");

            if (knownVariants == null || knownVariants.isEmpty()) {
                throw new UserException.CommandLineException(NO_DBSNP_EXCEPTION);
            }

            String referenceURL = referenceArguments.getReferenceFileName();

            if (readArguments.getReadFilesNames().size()!=1) {
                throw new UserException("Sorry, we only support a single reads input for now.");
            }
            String bam = readArguments.getReadFilesNames().get(0);

            if (!BucketUtils.isCloudStorageUrl(bam) && isRemote()) {
                throw new UserException("Sorry, for remote execution the BAM must be stored remotely.");
            }

            // Load the input bam
            final ReadsDataflowSource readsDataflowSource = new ReadsDataflowSource(bam, pipeline);
            readsHeader = readsDataflowSource.getHeader();
            final SAMSequenceDictionary readsDictionary = readsHeader.getSequenceDictionary();
            final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                                                       : IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());
            final PCollectionView<SAMFileHeader> headerSingleton = ReadsDataflowSource.getHeaderView(pipeline, readsHeader);
            final CountingReadFilter readFilterToApply = BaseRecalibrator.getStandardBQSRReadFilter(readsHeader);

            bunny.stepEnd("set up bam input");

            // Load the Variants and the Reference
            final ReferenceDataflowSource referenceDataflowSource = new ReferenceDataflowSource(pipeline.getOptions(), referenceURL, BQSR_REFERENCE_WINDOW_FUNCTION);
            bunny.stepEnd("create referenceDataflowSource");
            final SAMSequenceDictionary refDictionary = referenceDataflowSource.getReferenceSequenceDictionary(readsDictionary);
            bunny.stepEnd("load ref sequence dictionary");

            checkSequenceDictionaries(refDictionary, readsDictionary);
            bunny.stepEnd("checkSequenceDictionaries");
            PCollectionView<SAMSequenceDictionary> refDictionaryView = pipeline.apply(Create.of(refDictionary)).setName("refDictionary").apply(View.asSingleton());
            bunny.stepEnd("create ref dictionary view");


            List<SimpleInterval> shardedIntervals = IntervalUtils.cutToShards(intervals, bigShardSize);
            // since we currently can only read variants at the client, that's the right place to populate the shards.
            List<Variant> variants = VariantsDataflowSource.getVariantsList(knownVariants);

            bunny.stepEnd("load variants");
            ArrayList<AddContextDataToReadOptimized.ContextShard> shards = AddContextDataToReadOptimized.fillVariants(shardedIntervals, variants, margin);
            bunny.stepEnd("sharding variants");

            logger.info("Shipping "+shards.size()+" big shards.");

            PCollection<AddContextDataToReadOptimized.ContextShard> shardsPCol = pipeline.apply(Create.of(shards));

            PCollection<AddContextDataToReadOptimized.ContextShard> shardsWithContext = shardsPCol
                // big shards of variants -> smaller shards with variants, reads. We take the opportunity to filter the reads as close to the source as possible.
                .apply(ParDo.named("subdivideAndFillReads").of(AddContextDataToReadOptimized.subdivideAndFillReads(bam, smallShardSize, margin, readFilterToApply)))
                // add ref bases to the shards.
                .apply(ParDo.named("fillContext").of(AddContextDataToReadOptimized.fillContext(referenceDataflowSource)));

            final PCollection<RecalibrationTables> recalibrationTable = shardsWithContext.apply(new BaseRecalibratorOptimizedTransform(headerSingleton, refDictionaryView, recalArgs));

            if (null == serializedOutputTablesPath) {
                // we need those, so let's pick a temporary location for them.
                serializedOutputTablesPath = pickTemporaryRecaltablesPath(isRemote(), stagingLocation);
            }
            DataflowUtils.SaveDestination dest = DataflowUtils.serializeSingleObject(recalibrationTable, serializedOutputTablesPath);
            if (isRemote() && dest == DataflowUtils.SaveDestination.LOCAL_DISK) {
                throw new UserException("If running on the cloud, either leave serializedOutputTablesPath unset or point it to a GCS location.");
            }
            bunny.stepEnd("setup");
        } catch (UserException|GATKException rx) {
            throw rx;
        } catch (Exception x) {
            throw new GATKException("Unexpected: " + x.getMessage(), x);
        }
    }

    @Override
    protected void afterPipeline(Pipeline p) {
        bunny.stepEnd("dataflow");
        logger.info("Saving recalibration report to " + outputTablesPath);
        //  Get the table back and output it in text form to outputTablesPath
        // TODO: if running on the cloud and the output destination is on the cloud, then it's faster to have a worker do it directly, without the file roundtrip.
        try (ObjectInputStream oin = new ObjectInputStream(BucketUtils.openFile(serializedOutputTablesPath, p.getOptions()))) {
            Object o = oin.readObject();
            RecalibrationTables rt = (RecalibrationTables) o;
            BaseRecalibratorFn.saveTextualReport(new File(outputTablesPath), readsHeader, rt, recalArgs);
            bunny.stepEnd("repatriate_report");
        } catch (Exception e) {
            throw new GATKException("Unexpected: unable to read results file. (bug?)", e);
        }
        bunny.end();
    }

    protected static String pickTemporaryRecaltablesPath(boolean remote, String stagingLocation) {
        if (remote) {
            // in case multiple tests are running at the same time.
            String uniquifier = "-DF2-" + Math.abs(new Random().nextInt());
            return GcsPath.fromUri(stagingLocation).resolve(TEMP_RECALTABLES + uniquifier + ".sj").toString();
        } else {
            File outputTables = IOUtils.createTempFile(TEMP_RECALTABLES, ".sj");
            return outputTables.getPath();
        }
    }

    private void checkSequenceDictionaries(final SAMSequenceDictionary refDictionary, SAMSequenceDictionary readsDictionary) {
        Utils.nonNull(refDictionary);
        Utils.nonNull(readsDictionary);
        SequenceDictionaryUtils.validateDictionaries("reference", refDictionary, "reads", readsDictionary);
    }
}
