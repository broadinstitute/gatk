package org.broadinstitute.hellbender.tools.coveragemodel.germline;

import com.google.common.collect.ImmutableMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelArgumentCollection;
import org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelEMAlgorithm;
import org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelEMWorkspace;
import org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelParameters;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.IntegerCopyNumberTransitionMatrixCollection;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.IntegerCopyNumberTransitionProbabilityCacheCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.ContigGermlinePloidyAnnotationTableReader;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.GermlinePloidyAnnotatedTargetCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeDataCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeTableReader;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.*;
import java.util.ArrayList;
import java.util.Map;

/**
 * The command line tool for modelling and denoising coverage profiles and detecting germline copy number variation
 *
 * The tool can be run in two modes by specifying the job type (--jobType):
 *
 * <dl>
 *     <dt>LEARN_AND_CALL</dt>
 *     <dd>Learns model parameters from coverage profiles, saves the model, as well as various
 *     posteriors (copy number, sample-specific biases, read depth, model-compatibility log likelihoods, etc.)
 *     to disk</dd>
 *
 *     <dt>CALL_ONLY</dt>
 *     <dd>Takes a previously learned model from the argument --inputModelPath, and calculates various posteriors
 *     and saves the posteriors to disk</dd>
 * </dl>
 *
 * The tool requires the following mandatory arguments:
 *
 * <dl>
 *     <dt>Job type (specified via argument --jobType)</dt>
 *
 *     <dt>Combined raw or GC corrected (but not proportional) read counts table (specified via argument
 *     --input). The format is described in {@link ReadCountCollectionUtils}</dt>
 *
 *     <dt>Sample sex genotypes table (specified via argument --sexGenotypeTable). The format is described
 *     in {@link SexGenotypeTableReader}</dt>
 *
 *     <dt>Germline contig ploidy annotations for all sex genotypes (specified via argument
 *     --contigAnnotationsTable). The format is described in {@link ContigGermlinePloidyAnnotationTableReader}</dt>
 *
 *     <dt>Prior transition probabilities between different copy number states (specified via argument
 *     --copyNumberTransitionPriorTable). The format is described in
 *     {@link IntegerCopyNumberTransitionMatrixCollection.IntegerCopyNumberTransitionMatrixCollectionReader}</dt>
 *
 *     <dt>Output path for inferred model parameters, posteriors, and checkpoints (specified via argument
 *     --outputPath) </dt>
 * </dl>
 *
 * The tool will automatically use Spark clusters if available. Otherwise, it will run in the single-machine
 * (local) mode. If the user intends to run the tool on a single machine, it is recommended to disable Spark
 * altogether (--disableSpark true) since a local Spark context will only add unnecessary overhead.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Performs coverage profile modelling, denoising, and detecting germline copy number variation",
        oneLineSummary = "Performs coverage profile modelling, denoising, and detecting germline copy number variation",
        programGroup = CopyNumberProgramGroup.class
)
public final class GermlineCNVCaller extends SparkToggleCommandLineProgram {

    private static final long serialVersionUID = -1149969750485027900L;

    public enum JobType {
        /**
         * Learn model parameters and calculate posteriors
         */
        LEARN_AND_CALL,

        /**
         * Take a previously-learned set of model parameters and calculate posteriors
         */
        CALL_ONLY
    }

    private static final Logger logger = LogManager.getLogger(GermlineCNVCaller.class);

    public static final String FINAL_MODEL_SUBDIR = "model_final";
    public static final String FINAL_POSTERIORS_SUBDIR = "posteriors_final";

    public static final String CONTIG_PLOIDY_ANNOTATIONS_TABLE_LONG_NAME = "contigAnnotationsTable";
    public static final String CONTIG_PLOIDY_ANNOTATIONS_TABLE_SHORT_NAME = "CAT";

    public static final String SAMPLE_SEX_GENOTYPE_TABLE_LONG_NAME = "sexGenotypeTable";
    public static final String SAMPLE_SEX_GENOTYPE_TABLE_SHORT_NAME = "SGT";

    public static final String COPY_NUMBER_TRANSITION_PRIOR_TABLE_LONG_NAME = "copyNumberTransitionPriorTable";
    public static final String COPY_NUMBER_TRANSITION_PRIOR_TABLE_SHORT_NAME = "CNT";

    public static final String OUTPUT_PATH_LONG_NAME = "outputPath";
    public static final String OUTPUT_PATH_SHORT_NAME = "O";

    public static final String INPUT_MODEL_PATH_LONG_NAME = "inputModelPath";
    public static final String INPUT_MODEL_PATH_SHORT_NAME = "IMP";

    public static final String JOB_TYPE_LONG_NAME = "jobType";
    public static final String JOB_TYPE_SHORT_NAME = "JT";

    public static final String INPUT_READ_COUNTS_TABLE_LONG_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String INPUT_READ_COUNTS_TABLE_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    @Argument(
            doc = "Combined read count collection URI",
            fullName = INPUT_READ_COUNTS_TABLE_LONG_NAME,
            shortName = INPUT_READ_COUNTS_TABLE_SHORT_NAME,
            optional = false
    )
    protected String readCountsURI;

    @Argument(
            doc = "Contig ploidy annotations URI",
            fullName = CONTIG_PLOIDY_ANNOTATIONS_TABLE_LONG_NAME,
            shortName = CONTIG_PLOIDY_ANNOTATIONS_TABLE_SHORT_NAME,
            optional = false
    )
    protected String contigPloidyAnnotationsURI;

    @Argument(
            doc = "Sample sex genotypes URI",
            fullName = SAMPLE_SEX_GENOTYPE_TABLE_LONG_NAME,
            shortName = SAMPLE_SEX_GENOTYPE_TABLE_SHORT_NAME,
            optional = false
    )
    protected String sampleSexGenotypesURI;

    @Argument(
            doc = "Copy number transition prior table URI",
            fullName = COPY_NUMBER_TRANSITION_PRIOR_TABLE_LONG_NAME,
            shortName = COPY_NUMBER_TRANSITION_PRIOR_TABLE_SHORT_NAME,
            optional = false
    )
    protected String copyNumberTransitionPriorTableURI;

    @Argument(
            doc = "Output path for saving the model, posteriors, and checkpoints",
            fullName = OUTPUT_PATH_LONG_NAME,
            shortName = OUTPUT_PATH_SHORT_NAME,
            optional = false
    )
    protected String outputPath;

    @Argument(
            doc = "Job type",
            fullName = JOB_TYPE_LONG_NAME,
            shortName = JOB_TYPE_SHORT_NAME,
            optional = false
    )
    protected JobType jobType;

    @Argument(
            doc = "Input model (panel of normals) path",
            fullName = INPUT_MODEL_PATH_LONG_NAME,
            shortName = INPUT_MODEL_PATH_SHORT_NAME,
            optional = true
    )
    protected String modelPath = null;

    @ArgumentCollection
    protected final CoverageModelArgumentCollection params = new CoverageModelArgumentCollection();

    @ArgumentCollection
    protected final TargetArgumentCollection optionalTargets = new TargetArgumentCollection();

    /* custom serializers */
    private static final Map<String, String> coverageModellerExtraSparkProperties = ImmutableMap.of(
            "spark.kryo.registrator",
            "org.nd4j.Nd4jRegistrator,org.broadinstitute.hellbender.utils.spark.UnmodifiableCollectionsRegistrator");

    /**
     * Override doWork to inject custom Nd4j serializer and set a temporary checkpointing path
     */
    @Override
    protected Object doWork() {
        /* validate parameters */
        params.validate();

        JavaSparkContext ctx = null;
        if (!isDisableSpark) {
            /* create the spark context */
            final Map<String, String> sparkProperties = sparkArgs.getSparkProperties();
            appendExtraSparkProperties(sparkProperties, coverageModellerExtraSparkProperties);
            ctx = SparkContextFactory.getSparkContext(getProgramName(), sparkProperties, sparkArgs.getSparkMaster());
            ctx.setCheckpointDir(params.getRDDCheckpointingPath());
        } else {
            logger.info("Spark disabled.  sparkMaster option (" + sparkArgs.getSparkMaster() + ") ignored.");
        }

        try {
            runPipeline(ctx);
            return null;
        } finally {
            afterPipeline(ctx);
        }
    }

    /**
     * Checks {@param originalProperties} for keys present in {@param extraProperties}. For new keys,
     * it adds the value. For existing keys, it appends the value in a comma-separated fashion.
     *
     * @param originalProperties a non-null key-value {@link Map}
     * @param extraProperties a non-null key-value {@link Map}
     */
    private void appendExtraSparkProperties(@Nonnull final Map<String, String> originalProperties,
                                            @Nonnull final Map<String, String> extraProperties) {
        extraProperties.keySet().forEach(key ->
                originalProperties.put(key, originalProperties.containsKey(key)
                        ? originalProperties.get(key) + "," + extraProperties.get(key)
                        : extraProperties.get(key)));
    }

    /**
     * The main routine
     *
     * @param ctx a nullable Spark context
     */
    @Override
    protected void runPipeline(@Nullable JavaSparkContext ctx) {
        final TargetCollection<Target> optionalTargetsCollections = optionalTargets.readTargetCollection(true);
        if (optionalTargetsCollections == null) {
            logger.info("No target file was provided: using all targets in the combined read counts table");
        }

        logger.info("Parsing the read counts table...");
        final ReadCountCollection readCounts = loadReadCountCollection(optionalTargetsCollections);

        logger.info("Parsing the sample sex genotypes table...");
        final SexGenotypeDataCollection sexGenotypeDataCollection = loadSexGenotypeDataCollection();

        logger.info("Parsing the germline contig ploidy annotation table...");
        final GermlinePloidyAnnotatedTargetCollection ploidyAnnotatedTargetCollection =
                loadGermlinePloidyAnnotatedTargetCollection(readCounts);

        logger.info("Parsing the copy number transition prior table and initializing the caches...");
        final IntegerCopyNumberTransitionProbabilityCacheCollection transitionProbabilityCacheCollection =
                createIntegerCopyNumberTransitionProbabilityCacheCollection();

        final IntegerCopyNumberExpectationsCalculator integerCopyNumberExpectationsCalculator =
                new IntegerCopyNumberExpectationsCalculator(transitionProbabilityCacheCollection,
                        params.getMinLearningReadCount());

        final CoverageModelParameters model = getCoverageModelParameters();
        Utils.validateArg(model != null || !jobType.equals(JobType.CALL_ONLY),
                "Model parameters are not given; can not run the tool in the CALL_ONLY mode.");

        logger.info("Initializing the EM algorithm workspace...");
        final IntegerCopyNumberReferenceStateFactory referenceStateFactory =
                new IntegerCopyNumberReferenceStateFactory(ploidyAnnotatedTargetCollection);

        final CoverageModelEMWorkspace<IntegerCopyNumberState> workspace = new CoverageModelEMWorkspace<>(
                readCounts, ploidyAnnotatedTargetCollection, sexGenotypeDataCollection,
                integerCopyNumberExpectationsCalculator, params, model, referenceStateFactory, ctx);
        final CoverageModelEMAlgorithm<IntegerCopyNumberState> algo = new CoverageModelEMAlgorithm<>(params,
                workspace);

        switch (jobType) {
            case LEARN_AND_CALL:
                algo.runExpectationMaximization();
                logger.info("Saving the model to disk...");
                workspace.writeModel(new File(outputPath, FINAL_MODEL_SUBDIR).getAbsolutePath());
                break;

            case CALL_ONLY:
                algo.runExpectation();
                break;

            default:
                throw new UnsupportedOperationException(String.format("\"%s\" is not recognized as a supported job type",
                        jobType.name()));
        }

        logger.info("Saving posteriors to disk...");
        workspace.writePosteriors(new File(outputPath, FINAL_POSTERIORS_SUBDIR).getAbsolutePath(),
                CoverageModelEMWorkspace.PosteriorVerbosityLevel.EXTENDED);
    }

    private CoverageModelParameters getCoverageModelParameters() {
        CoverageModelParameters model;
        if (modelPath != null) {
            logger.info("Loading model parameters...");
            model = CoverageModelParameters.read(modelPath);
        } else {
            model = null;
        }
        return model;
    }

    private IntegerCopyNumberTransitionProbabilityCacheCollection createIntegerCopyNumberTransitionProbabilityCacheCollection() {
        IntegerCopyNumberTransitionProbabilityCacheCollection transitionProbabilityCacheCollection;
        try (final Reader copyNumberTransitionPriorTableReader = getReaderFromURI(copyNumberTransitionPriorTableURI)) {
            final String parentURI = new File(copyNumberTransitionPriorTableURI).getParent();
            final IntegerCopyNumberTransitionMatrixCollection transitionMatrixCollection =
                    IntegerCopyNumberTransitionMatrixCollection.read(copyNumberTransitionPriorTableReader, parentURI);
            transitionProbabilityCacheCollection = new IntegerCopyNumberTransitionProbabilityCacheCollection(
                    transitionMatrixCollection, true);
        } catch (final IOException ex) {
            throw new RuntimeException("Could not close the copy number transition prior table reader", ex);
        }
        return transitionProbabilityCacheCollection;
    }

    private GermlinePloidyAnnotatedTargetCollection loadGermlinePloidyAnnotatedTargetCollection(@Nonnull final ReadCountCollection readCounts) {
        GermlinePloidyAnnotatedTargetCollection ploidyAnnotatedTargetCollection;
        try (final Reader ploidyAnnotationsReader = getReaderFromURI(contigPloidyAnnotationsURI)) {
            ploidyAnnotatedTargetCollection = new GermlinePloidyAnnotatedTargetCollection(ContigGermlinePloidyAnnotationTableReader
                    .readContigGermlinePloidyAnnotationsFromReader(contigPloidyAnnotationsURI, ploidyAnnotationsReader),
                    new ArrayList<>(readCounts.targets()));
        } catch (final IOException ex) {
            throw new RuntimeException("Could not close the germline contig ploidy annotations table reader", ex);
        }
        return ploidyAnnotatedTargetCollection;
    }

    private SexGenotypeDataCollection loadSexGenotypeDataCollection() {
        SexGenotypeDataCollection sexGenotypeDataCollection;
        try (final Reader sexGenotypeDataCollectionReader = getReaderFromURI(sampleSexGenotypesURI)) {
            sexGenotypeDataCollection = new SexGenotypeDataCollection(sexGenotypeDataCollectionReader,
                    sampleSexGenotypesURI);
        } catch (final IOException ex) {
            throw new RuntimeException("Could not close the input sample sex genotypes table reader", ex);
        }
        return sexGenotypeDataCollection;
    }

    private ReadCountCollection loadReadCountCollection(@Nullable final TargetCollection<Target> targetsCollections) {
        ReadCountCollection readCounts;
        try (final Reader readCountsReader = getReaderFromURI(readCountsURI)) {
            if (targetsCollections == null) {
                readCounts = ReadCountCollectionUtils.parse(readCountsReader, readCountsURI);
            } else {
                readCounts = ReadCountCollectionUtils.parse(readCountsReader, readCountsURI, targetsCollections, true);
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not parse the read counts table", ex);
        }
        return readCounts;
    }

    /**
     * Takes a URI string (a local file, a Google bucket file, or an HDFS file) and returns a {@link Reader}
     *
     * @param path input URI string
     * @return an instance of {@link Reader}
     */
    private Reader getReaderFromURI(@Nonnull final String path) {
        final InputStream inputStream = BucketUtils.openFile(path);
        return new BufferedReader(new InputStreamReader(inputStream));
    }
}
