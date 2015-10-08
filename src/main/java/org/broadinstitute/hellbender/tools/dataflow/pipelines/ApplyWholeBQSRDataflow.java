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
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ContextShard;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsShard;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantsDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.transforms.composite.AddContextDataToRead;
import org.broadinstitute.hellbender.engine.dataflow.transforms.composite.AddContextDataToReadOptimized;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.ApplyBQSRTransformOptimized;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalOutput;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorFn;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorOptimizedTransform;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.logging.BunnyLog;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.File;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


/**
 * Optimized BaseRecalibratorDataflow + ApplyBQSR.
 * I expect this to go away once we fix markdups and can move this code into the pipeline
 */
@CommandLineProgramProperties(
        summary = "BaseRecalibrator + ApplyBQSR",
        oneLineSummary = "Generates recalibration table",
        programGroup = ReadProgramGroup.class
)
public class ApplyWholeBQSRDataflow extends DataflowCommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
      * Reference window function for BQSR. For each read, returns an interval representing the span of
      * reference bases required by the BQSR algorithm for that read. Should be passed into the
      * {@link ReferenceDataflowSource} object for the {@link AddContextDataToRead} transform.
      */
    public static final SerializableFunction<GATKRead, SimpleInterval> BQSR_REFERENCE_WINDOW_FUNCTION =
            read -> BAQ.getReferenceWindowForRead(read, BAQ.DEFAULT_BANDWIDTH, BAQ.DEFAULT_INCLUDE_CLIPPED_BASES);

    public static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";

    private final static Logger logger = LogManager.getLogger(ApplyWholeBQSRDataflow.class);
    // temporary file with the serialized recalibrationTables.
    private final static String TEMP_RECALTABLES = "temp-ds-recaltables";

    // ------------------------------------------
    // Command-line options

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private transient final RecalibrationArgumentCollection BRAC = new RecalibrationArgumentCollection();

    @ArgumentCollection
    private final RequiredReadInputArgumentCollection readArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    private final IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @ArgumentCollection
    private final ReferenceInputArgumentCollection referenceArguments = new RequiredReferenceInputArgumentCollection();

    @Argument(doc = "the known variants", shortName = "knownSites", fullName = "knownSites", optional = false)
    private List<String> baseRecalibrationKnownVariants;

    /**
     * Path to save the final recalibration tables to. Local.
     */
    @Argument(doc = "Path to save the recalibrationTables to.",
            fullName = "RECAL_TABLE_FILE", optional = true)
    protected String outputTablesPath = null;

    /**
     * Path to save the temporary serialized recalibration tables to. Local or GCS.
     */
    @Argument(doc = "Path to save the serialized recalibrationTables to. If running on the cloud, either leave serializedOutputTablesPath unset or point it to a GCS location.",
            shortName = "sotp", fullName = "serializedOutputTablesPath", optional = true)
    protected String serializedOutputTablesPath = null;

    private transient final ApplyBQSRArgumentCollection ABAC = new ApplyBQSRArgumentCollection();

    @Argument(fullName="largeShard", optional=true)
    public int largeShard = 1_000_000;

    @Argument(fullName="smallShard", optional=true)
    public int smallShard = 5_000;

    @Argument(fullName="jobNameExtra", optional=true)
    public String jobNameExtra = null;
    // ------------------------------------------

    // Whether we want to save the textual version of the output.
    private boolean saveTextualTables;

    private SAMFileHeader readsHeader;

    private transient BunnyLog bunny = new BunnyLog(logger);

    @Override
    protected String jobName() {
        String s = "";
        if (jobNameExtra!=null) s = jobNameExtra+"-";
        return "ApplyWholeBQSR-"+ s + (100+new Random().nextInt(900));
    }

    @Override
    protected void setupPipeline(Pipeline pipeline) {
        try {
            bunny.start("ApplyWholeBQSR");

            if (null==baseRecalibrationKnownVariants || baseRecalibrationKnownVariants.size()==0) {
                throw new UserException.CommandLineException(NO_DBSNP_EXCEPTION);
            }

            // no saving for now.
            saveTextualTables = false;

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
            CountingReadFilter filterToApply = ApplyWholeBQSRDataflow.readFilter(readsHeader);

            bunny.stepEnd("set up bam input");

            // Load the Variants and the Reference
            //final VariantsDataflowSource variantsDataflowSource = new VariantsDataflowSource(baseRecalibrationKnownVariants, pipeline);
            final ReferenceDataflowSource referenceDataflowSource = new ReferenceDataflowSource(pipeline.getOptions(), referenceURL, BQSR_REFERENCE_WINDOW_FUNCTION);
            bunny.stepEnd("create referenceDataflowSource");
            final SAMSequenceDictionary refDictionary = referenceDataflowSource.getReferenceSequenceDictionary(readsDictionary);
            bunny.stepEnd("load ref sequence dictionary");

            checkSequenceDictionaries(refDictionary, readsDictionary);
            bunny.stepEnd("checkSequenceDictionaries");
            PCollectionView<SAMSequenceDictionary> refDictionaryView = pipeline.apply(Create.of(refDictionary)).setName("refDictionary").apply(View.asSingleton());
            bunny.stepEnd("create ref dictionary view");

            // the granularity at which we'll want to assign work and read inputs.
            final int bigShardSize = largeShard;
            // the granularity at which we want to batch reads,variants,and reference for processing.
            final int smallShardSize = smallShard;

            List<SimpleInterval> shardedIntervals = IntervalUtils.cutToShards(intervals, bigShardSize);
            // since we currently can only read variants at the client, that's the right place to populate the shards.
            List<Variant> variants = VariantsDataflowSource.getVariantsList(baseRecalibrationKnownVariants);

            bunny.stepEnd("load variants");
            int margin = 1000;
            ArrayList<ContextShard> shards = fillVariants(shardedIntervals, variants, margin);
            bunny.stepEnd("sharding variants");

            logger.info("Shipping "+shards.size()+" big shards.");

            PCollection<ContextShard> shardsPCol = pipeline.apply(Create.of(shards)).setName("variants");

            PCollection<ContextShard> shardsNoContext = shardsPCol
                // big shards of variants -> smaller shards with variants, reads. We take the opportunity to filter the reads as close to the source as possible.
                .apply(ParDo.named("subdivideAndFillReads").of(AddContextDataToReadOptimized.subdivideAndFillReads(bam, smallShardSize, margin, null)));
            PCollection<ReadsShard> shardsReadOnly = shardsNoContext
                .apply(ParDo.named("extract reads").of(AddContextDataToReadOptimized.extractReads()));
            PCollection<ContextShard> shardsWithContext = shardsNoContext
                // add ref bases to the shards
                .apply(ParDo.named("fillContext").of(AddContextDataToReadOptimized.fillContext(referenceDataflowSource, filterToApply))
                        .withSideInputs(refDictionaryView));

            BaseRecalibratorOptimizedTransform baseRecal = new BaseRecalibratorOptimizedTransform(headerSingleton, refDictionaryView, BRAC);
            final PCollection<BaseRecalOutput> baseRecalOutput = shardsWithContext
                .apply(baseRecal)
                .apply(baseRecal.toBaseRecalOutput());

            final PCollectionView<BaseRecalOutput> baseRecalOutputView = baseRecalOutput.apply(View.asSingleton());

            bunny.stepEnd("baseRecal setup");

            //final PCollectionView<BaseRecalOutput> recalInfoSingletonView = BaseRecalOutputSource.loadFileOrRemote(pipeline, BQSR_RECAL_FILE_NAME).apply(View.asSingleton());

            final PCollection<GATKRead> output =  shardsReadOnly
                .apply(new ApplyBQSRTransformOptimized(headerSingleton, baseRecalOutputView, ABAC));

            // uncomment to save the recalibration tables for debug
            /*
            // If saving textual output then we need to make sure we can get to the output
            if (saveTextualTables) {
                if (null == serializedOutputTablesPath) {
                    // we need those, so let's pick a temporary location for them.
                    outputTablesPath = pickTemporaryRecaltablesPath(isRemote(), stagingLocation);
                }
                DataflowUtils.SaveDestination dest = DataflowUtils.serializeSingleObject(recalibrationTable, outputTablesPath);
                if (isRemote() && dest == DataflowUtils.SaveDestination.LOCAL_DISK) {
                    throw new UserException("If running on the cloud, either leave outputTablesPath unset or point it to a GCS location.");
                }
            }
            */
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

        if (saveTextualTables) {
            logger.info("Saving recalibration report to "+outputTablesPath);
            //  Get the table back and output it in text form to RAC.RECAL_TABLE.
            // TODO: if running on the cloud and the output destination is on the cloud, then it's faster to have a worker do it directly, without the file roundtrip.
            try (ObjectInputStream oin = new ObjectInputStream(BucketUtils.openFile(serializedOutputTablesPath, p.getOptions()))) {
                Object o = oin.readObject();
                RecalibrationTables rt = (RecalibrationTables) o;
                BaseRecalibratorFn.saveTextualReport(new File(outputTablesPath), readsHeader, rt, BRAC);
                bunny.stepEnd("repatriate_report");
            } catch (Exception e) {
                throw new GATKException("Unexpected: unable to read results file. (bug?)", e);
            }
        } else {
            // because the user didn't ask us to.
            logger.info("Not saving recalibration report");
        }
        bunny.end();
    }

    protected static String pickTemporaryRecaltablesPath(boolean remote, String stagingLocation) {
        if (remote) {
            return GcsPath.fromUri(stagingLocation).resolve(TEMP_RECALTABLES + ".sj").toString();
        } else {
            File outputTables = IOUtils.createTempFile(TEMP_RECALTABLES, ".sj");
            return outputTables.getPath();
        }
    }

    private static CountingReadFilter readFilter( final SAMFileHeader header ) {
        return new CountingReadFilter("Wellformed", new WellformedReadFilter(header))
            .and(new CountingReadFilter("Mapping_Quality_Not_Zero", ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO))
            .and(new CountingReadFilter("Mapping_Quality_Available", ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE))
            .and(new CountingReadFilter("Mapped", ReadFilterLibrary.MAPPED))
            .and(new CountingReadFilter("Primary_Alignment", ReadFilterLibrary.PRIMARY_ALIGNMENT))
            .and(new CountingReadFilter("Not_Duplicate", ReadFilterLibrary.NOT_DUPLICATE))
            .and(new CountingReadFilter("Passes_Vendor_Quality_Check", ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK));
    }

    private void checkSequenceDictionaries(final SAMSequenceDictionary refDictionary, SAMSequenceDictionary readsDictionary) {
        Utils.nonNull(refDictionary);
        Utils.nonNull(readsDictionary);
        SequenceDictionaryUtils.validateDictionaries("reference", refDictionary, "reads", readsDictionary);
    }

    // Given a list of shards and a list of variants,
    // add each variant to every (shard+margin) that it overlaps
    private ArrayList<ContextShard> fillVariants(List<SimpleInterval> shardedIntervals, List<Variant> variants, int margin) {
        ArrayList<ContextShard> ret = new ArrayList<>();
        IntervalsSkipList<Variant> fastVars = new IntervalsSkipList<>(variants);
        for (SimpleInterval s : shardedIntervals) {
            int start = s.getStart() - margin;
            if (start<1) start=1;
            int end = s.getEnd() + margin;
            // here it's OK if end is past the contig's boundary, there just won't be any variant there.
            SimpleInterval interval = new SimpleInterval(s.getContig(), start, end);
            ret.add(new ContextShard(s).withVariants(fastVars.getOverlapping(interval)));
        }
        return ret;
    }

}
