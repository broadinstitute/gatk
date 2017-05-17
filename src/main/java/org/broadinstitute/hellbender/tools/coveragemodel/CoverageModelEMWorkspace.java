package org.broadinstitute.hellbender.tools.coveragemodel;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.AbstractUnivariateSolver;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.storage.StorageLevel;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelCopyRatioEmissionProbabilityCalculator.EmissionCalculationStrategy;
import org.broadinstitute.hellbender.tools.coveragemodel.interfaces.CopyRatioExpectationsCalculator;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.FourierLinearOperatorNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.GeneralLinearOperator;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.IterativeLinearSolverNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.IterativeLinearSolverNDArray.ExitStatus;
import org.broadinstitute.hellbender.tools.coveragemodel.math.RobustBrentSolver;
import org.broadinstitute.hellbender.tools.coveragemodel.math.SynchronizedUnivariateSolver;
import org.broadinstitute.hellbender.tools.coveragemodel.math.UnivariateSolverSpecifications;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.GermlinePloidyAnnotatedTargetCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeData;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeDataCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.interfaces.AlleleMetadataProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.CallStringProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.ScalarProducer;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HMMSegmentProcessor;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecord;
import org.broadinstitute.hellbender.utils.hmm.segmentation.HiddenStateSegmentRecordWriter;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.buffer.DataBuffer;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.INDArrayIndex;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.nd4j.linalg.ops.transforms.Transforms;
import scala.Tuple2;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * This class represents the driver-node workspace for EM algorithm calculations of the coverage model.
 *
 * At the moment, this class is responsible for:
 *
 *   - initializing the model parameters
 *   - initializing workers
 *   - performing various E-step and M-step updates (called by {@link CoverageModelEMAlgorithm}
 *   - writing posteriors to disk
 *
 * The user instantiates this class and passes it on to {@link CoverageModelEMAlgorithm}; the latter calls
 * various E-step and M-step methods. Once the EM algorithm is converged, the user can write out the latest
 * posteriors to disk by calling {@link #writePosteriors}.
 *
 * TODO github/gatk-protected issue #1021 -- this class must be refactored; much of the exposed public
 * methods must be made private
 *
 * @param <STATE> copy ratio (or copy number) hidden state type
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelEMWorkspace<STATE extends AlleleMetadataProducer & CallStringProducer &
        ScalarProducer> {

    private static final Logger logger = LogManager.getLogger(CoverageModelEMWorkspace.class);

    protected final CoverageModelArgumentCollection config;

    /**
     * The input read count collection after initial processing
     */
    protected final ReadCountCollection processedReadCounts;

    /**
     * Target list after initial processing
     */
    protected final List<Target> processedTargetList;

    /**
     * Sample names list after initial processing
     */
    protected final List<String> processedSampleNameList;

    /**
     * Sex genotype data list after initial processing
     */
    protected final List<SexGenotypeData> processedSampleSexGenotypeData;

    /**
     * Model parameters adapted to the targets in {@link CoverageModelEMWorkspace#processedReadCounts}
     */
    protected final CoverageModelParameters processedModel;

    /**
     * An implementation of {@link CopyRatioExpectationsCalculator} for calculating copy ratio (or copy number)
     * prior and posterior expectations
     */
    protected final CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, STATE> copyRatioExpectationsCalculator;

    /**
     * Collection of sex genotype data of the input samples
     */
    protected final SexGenotypeDataCollection sexGenotypeDataCollection;

    /**
     * Collection of germline ploidies of the targets for each sex genotype
     */
    protected final GermlinePloidyAnnotatedTargetCollection germlinePloidyAnnotatedTargetCollection;

    /**
     * Number of samples in {@link CoverageModelEMWorkspace#processedReadCounts}
     */
    protected final int numSamples;

    /**
     * Number of targets in {@link CoverageModelEMWorkspace#processedReadCounts}
     */
    protected final int numTargets;

    /**
     * Dimension of the bias latent space
     */
    protected final int numLatents;

    /**
     * Whether or not a Spark context is available
     */
    protected final boolean sparkContextIsAvailable;

    /**
     * The Spark context (null in the local mode)
     */
    protected final JavaSparkContext ctx;

    private final BiFunction<SexGenotypeData, Target, STATE> referenceStateFactory;

    /* Spark-related members --- BEGIN */

    /**
     * The latest RDD of compute blocks
     */
    private JavaPairRDD<LinearlySpacedIndexBlock, CoverageModelEMComputeBlock> computeRDD;

    /**
     * The latest checkpointed RDD of compute blocks
     */
    private JavaPairRDD<LinearlySpacedIndexBlock, CoverageModelEMComputeBlock> prevCheckpointedComputeRDD;

    /**
     * A deque of cached RDDs of compute blocks
     */
    private Deque<JavaPairRDD<LinearlySpacedIndexBlock, CoverageModelEMComputeBlock>> prevCachedComputeRDDDeque = new LinkedList<>();

    /**
     * Counts the number of calls made to {@link CoverageModelEMWorkspace#cacheWorkers(String)} since
     * the most recent RDD checkpointing call
     */
    private int cacheCallCounter;

    /* END --- Spark-related members */

    /**
     * An instance of {@link CoverageModelEMComputeBlock} for the local mode
     */
    private CoverageModelEMComputeBlock localComputeBlock;

    /**
     * Number of target-space blocks
     */
    protected final int numTargetBlocks;

    /**
     * List of target-space block specs
     */
    protected List<LinearlySpacedIndexBlock> targetBlocks;

    /* Driver node copy of posteriors --- BEGIN */

    /**
     * $E[log(d_s)]$
     */
    protected final INDArray sampleMeanLogReadDepths;

    /**
     * $Var[log(d_s)]$
     */
    protected final INDArray sampleVarLogReadDepths;

    /**
     * $E[z_{sm}]$
     */
    protected final INDArray sampleBiasLatentPosteriorFirstMoments;

    /**
     * $E[z_{sm} z_{sn}]$
     */
    protected final INDArray sampleBiasLatentPosteriorSecondMoments;

    /**
     * $E[\gamma_s]$
     */
    protected final INDArray sampleUnexplainedVariance;

    /**
     * Posterior expectation of HMM backbone (log likelihood - log emission likelihood)
     */
    protected final INDArray sampleLogChainPosteriors;

    /**
     * Norm_2 of bias covariates
     */
    protected final INDArray biasCovariatesNorm2;

    /* END --- Driver node copy of posteriors */

    /* Driver node copy of model parameters --- BEGIN */

    /**
     * Bias covariates ARD coefficients
     */
    protected final INDArray biasCovariatesARDCoefficients;

    protected final List<INDArray> biasCovariatesARDCoefficientsHistory;

    /* END --- Driver node copy of model parameters */

    /**
     * Fourier factors of the CNV-avoiding regularizer
     */
    private final double[] fourierFactors;

    protected final boolean biasCovariatesEnabled;

    private final boolean ardEnabled;

    protected final List<Double> logLikelihoodHistory;

    /**
     * Annotation for methods that call cache() on an RDD
     */
    private @interface CachesRDD {
    }

    /**
     * Annotation for methods that evaluate an RDD
     */
    private @interface EvaluatesRDD {
    }

    /**
     * Annotation for methods that map an RDD
     */
    private @interface UpdatesRDD {
    }

    /**
     * Public constructor
     *
     * @param rawReadCounts an instance of {@link ReadCountCollection} containing raw read counts
     * @param germlinePloidyAnnotatedTargetCollection an instance of {@link GermlinePloidyAnnotatedTargetCollection}
     *                                                for obtaining target ploidies for different sex genotypes
     * @param sexGenotypeDataCollection an instance of {@link SexGenotypeDataCollection} for obtaining sample sex genotypes
     * @param config coverage model EM algorithm configuration
     * @param ctx the Spark context
     * @param copyRatioExpectationsCalculator an implementation of {@link CopyRatioExpectationsCalculator}
     */
    @UpdatesRDD @CachesRDD @EvaluatesRDD
    public CoverageModelEMWorkspace(@Nonnull final ReadCountCollection rawReadCounts,
                                    @Nonnull final GermlinePloidyAnnotatedTargetCollection germlinePloidyAnnotatedTargetCollection,
                                    @Nonnull final SexGenotypeDataCollection sexGenotypeDataCollection,
                                    @Nonnull final CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, STATE> copyRatioExpectationsCalculator,
                                    @Nonnull final CoverageModelArgumentCollection config,
                                    @Nullable final CoverageModelParameters model,
                                    @Nonnull final BiFunction<SexGenotypeData, Target, STATE> referenceStateFactory,
                                    @Nullable final JavaSparkContext ctx) {
        this.config = Utils.nonNull(config, "Coverage model EM-algorithm parameters must be non-null");
        this.copyRatioExpectationsCalculator = Utils.nonNull(copyRatioExpectationsCalculator, "Copy ratio posterior calculator" +
                " must be non-null");
        this.germlinePloidyAnnotatedTargetCollection = Utils.nonNull(germlinePloidyAnnotatedTargetCollection,
                "The germline ploidy-annotated target collection must be non-null");
        this.sexGenotypeDataCollection = Utils.nonNull(sexGenotypeDataCollection,
                "The sex genotype data collection must be non-null");
        this.referenceStateFactory = Utils.nonNull(referenceStateFactory,
                "The reference state factory must be non-null");
        Utils.nonNull(rawReadCounts, "Raw read count collection must be non-null");

        /* perform basic consistency check of arguments */
        validateWorkspaceArgs(rawReadCounts, germlinePloidyAnnotatedTargetCollection, sexGenotypeDataCollection);

        /* check Nd4j data type */
        if (!Nd4j.dataType().equals(DataBuffer.Type.DOUBLE)) {
            throw new GATKException("Nd4j data type must be set to double for coverage modeller routines" +
                    " to function properly. This can be done by setting JVM system property \"dtype\" to" +
                    " \"double\". Can not continue.");
        }

        /* assert that targets are lexicographically sorted; otherwise, sort them */
        final ReadCountCollection targetSortedRawReadCounts = getTargetSortedReadCountCollection(rawReadCounts);

        /* if a model is provided, adapt the target list and number or latents to the model */
        if (model != null) {
            final ReadCountCollection intermediateReadCounts = processReadCountCollection(targetSortedRawReadCounts, config,
                    sexGenotypeDataCollection, germlinePloidyAnnotatedTargetCollection);
            /* adapt model and read count collection */
            final ImmutablePair<CoverageModelParameters, ReadCountCollection> modelReadCountsPair =
                    CoverageModelParameters.adaptModelToReadCountCollection(model, intermediateReadCounts, logger);
            processedModel = modelReadCountsPair.left;
            processedReadCounts = modelReadCountsPair.right;
            numLatents = processedModel.getNumLatents();
            Utils.validateArg(processedModel.getNumLatents() == config.getNumLatents(), "Discrepancy between the dimension" +
                    " of the bias latent space between the provided model and the arguments");
            Utils.validateArg(processedModel.isARDEnabled() == config.isARDEnabled(), "The model does not have ARD" +
                    " information -- ARD must be disabled in the tool arguments");
        } else {
            processedModel = null;
            processedReadCounts = processReadCountCollection(targetSortedRawReadCounts, config, sexGenotypeDataCollection,
                    germlinePloidyAnnotatedTargetCollection);
            numLatents = config.getNumLatents();
        }

        ardEnabled = config.isARDEnabled();
        biasCovariatesEnabled = numLatents > 0;

        numSamples = processedReadCounts.columnNames().size();
        numTargets = processedReadCounts.targets().size();
        processedTargetList = Collections.unmodifiableList(processedReadCounts.targets());
        processedSampleNameList = Collections.unmodifiableList(processedReadCounts.columnNames());
        processedSampleSexGenotypeData = Collections.unmodifiableList(processedSampleNameList.stream()
                .map(sexGenotypeDataCollection::getSampleSexGenotypeData)
                .collect(Collectors.toList()));
        logLikelihoodHistory = new ArrayList<>();

        logger.info(String.format("Number of samples before and after pre-processing read counts: (%d, %d)",
                rawReadCounts.columnNames().size(), numSamples));
        logger.info(String.format("Number of targets before and after pre-processing read counts: (%d, %d)",
                rawReadCounts.targets().size(), numTargets));

        this.ctx = ctx;
        if (ctx == null) {
            sparkContextIsAvailable = false;
            this.numTargetBlocks = 1;
        } else {
            sparkContextIsAvailable = true;
            this.numTargetBlocks = ParamUtils.inRange(config.getNumTargetSpacePartitions(), 1, numTargets,
                    "Number of target blocks must be between 1 and the size of target space.");
        }

        /* allocate memory and initialize driver-node copy of posteriors */
        sampleMeanLogReadDepths = Nd4j.zeros(numSamples, 1);
        sampleVarLogReadDepths = Nd4j.zeros(numSamples, 1);
        sampleUnexplainedVariance = Nd4j.zeros(numSamples, 1);
        sampleLogChainPosteriors = Nd4j.zeros(numSamples, 1);

        if (biasCovariatesEnabled) {
            sampleBiasLatentPosteriorFirstMoments = Nd4j.zeros(numSamples, numLatents);
            sampleBiasLatentPosteriorSecondMoments = Nd4j.zeros(numSamples, numLatents, numLatents);
            biasCovariatesNorm2 = Nd4j.zeros(1, numLatents);
        } else {
            logger.info("Bias covariates are disabled");
            sampleBiasLatentPosteriorFirstMoments = null;
            sampleBiasLatentPosteriorSecondMoments = null;
            biasCovariatesNorm2 = null;
        }

        /* allocate memory for driver-node copy of model parameters */
        if (ardEnabled) { /* if so, numLatents > 0 by way of pre-validating {@code config} */
            biasCovariatesARDCoefficients = Nd4j.zeros(new int[] {1, numLatents}).addi(config.getInitialARDPrecisionAbsolute());
            biasCovariatesARDCoefficientsHistory = new ArrayList<>();
        } else {
            logger.info("ARD for bias covariates is disabled");
            biasCovariatesARDCoefficients = null;
            biasCovariatesARDCoefficientsHistory = null;
        }

        /* initialize the CNV-avoiding regularizer filter */
        if (config.fourierRegularizationEnabled()) {
            fourierFactors = GATKProtectedMathUtils.getMidpassFilterFourierFactors(numTargets,
                    numTargets/config.getMaximumCNVLength(), numTargets/config.getMinimumCNVLength());
        } else {
            fourierFactors = null;
        }

        /* initialize target-space blocks */
        initializeTargetSpaceBlocks();

        /* create compute blocks */
        instantiateWorkers();

        /* push read counts to compute blocks and initialize copy ratio posterior expectations to prior expectations */
        pushInitialDataToComputeBlocks();

        /* initialize worker posteriors */
        initializeWorkerPosteriors();

        /* initialize model parameters */
        initializeWorkerModelParameters();
    }

    /**
     * Perform basic consistency check of arguments:
     *
     * <dl>
     *     <dt> target names in the read count collection is unique </dt>
     *     <dt> targets in the read count collection exists in the germline ploidy-annotated target collection </dt>
     *     <dt> all samples have sex genotype metadata </dt>
     *     <dt> contig germline ploidy annotations are available for all required sex genotyes </dt>
     * </dl>
     * @param rawReadCounts not {@code null} instance of {@link ReadCountCollection}
     * @param germlinePloidyAnnotatedTargetCollection an instance of {@link GermlinePloidyAnnotatedTargetCollection}
     *                                                for obtaining target ploidies for different sex genotypes
     * @param sexGenotypeDataCollection an instance of {@link SexGenotypeDataCollection} for obtaining sample sex genotypes
     */
    private void validateWorkspaceArgs(@Nonnull final ReadCountCollection rawReadCounts,
                                       @Nonnull final GermlinePloidyAnnotatedTargetCollection germlinePloidyAnnotatedTargetCollection,
                                       @Nonnull final SexGenotypeDataCollection sexGenotypeDataCollection) {
        /* assert that targets in the read count collection have unique names */
        Utils.validateArg(rawReadCounts.targets().stream().map(Target::getName).collect(Collectors.toSet()).size() ==
                rawReadCounts.targets().size(), "The targets in the read count collection must have unique names");

        /* assert that targets in the read count count collection have germline ploidy annotations */
        Utils.validateArg(rawReadCounts.targets().stream()
                        .allMatch(t -> germlinePloidyAnnotatedTargetCollection.getFullTargetSet().contains(t)),
                String.format("Some of the targets in the read count collection do not have germline contig ploidy"
                        + " annotations: %s", Sets.difference(rawReadCounts.targets().stream()
                                .map(Target::getName).collect(Collectors.toSet()),
                        germlinePloidyAnnotatedTargetCollection.getFullTargetSet().stream().map(Target::getName)
                                .collect(Collectors.toSet())).stream().collect(Collectors.joining(", "))));

        /* assert that all samples have sex genotype data */
        Utils.validateArg(rawReadCounts.columnNames().stream()
                        .allMatch(s -> sexGenotypeDataCollection.getSampleNames().contains(s)),
                String.format("Some of the samples in the read count collection do not have sex genotype annotations: %s",
                        Sets.difference(new HashSet<>(rawReadCounts.columnNames()),
                                sexGenotypeDataCollection.getSampleNames()).stream().collect(Collectors.joining(", "))));

        /* assert that all sex genotypes have contig ploidy annotations */
        final Set<String> sampleSexGenotypesSet = sexGenotypeDataCollection.getSexGenotypeDataList().stream()
                .map(SexGenotypeData::getSexGenotype).collect(Collectors.toSet());
        for (final Target t : rawReadCounts.targets()) {
            final Set<String> missingSexGenotypes = Sets.difference(
                    sampleSexGenotypesSet,
                    germlinePloidyAnnotatedTargetCollection.getContigGermlinePloidyAnnotation(t).getGenotypesSet());
            Utils.validateArg(missingSexGenotypes.isEmpty(), String.format("Some of the sample sex genotypes do not" +
                    " have germline contig ploidy annotations: %s", missingSexGenotypes.stream()
                    .collect(Collectors.joining(", "))));
        }
    }

    /**
     * If the targets in the read count collection are lexicographically sorted, does nothing (returns the original
     * reference). Otherwise, returns a copy of the read count collection with sorted targets.
     *
     * @param originalReadCountCollection a non-null instance of {@link ReadCountCollection}
     * @return an instance of {@link ReadCountCollection}
     */
    private static ReadCountCollection getTargetSortedReadCountCollection(@Nonnull final ReadCountCollection originalReadCountCollection) {
        final List<Target> originalTargetList = originalReadCountCollection.targets();
        if (!IntStream.range(0, originalTargetList.size()-1)
                .allMatch(ti -> IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR
                        .compare(originalTargetList.get(ti + 1), originalTargetList.get(ti)) > 0)) {
            final List<Target> sortedTargetList =  originalTargetList.stream()
                    .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());
            return originalReadCountCollection.arrangeTargets(sortedTargetList);
        } else {
            return originalReadCountCollection;
        }
    }

    /**
     * Process read counts and filter targets and/or samples as follows:
     *
     * <dl>
     *     <dt> Remove targets with bad values (NaN, infinity, negative) </dt>
     *     <dt> Remove targets that will be masked for learning across all samples </dt>
     * </dl>
     *
     * @implNote it is best to keep this filter static and functional in design for robustness/clarity.
     *
     * TODO github/gatk-protected issue #855 -- consider adding more filters
     * - remove targets with very high and very low GC content (can be done externally)
     * - remove targets with lots of repeats (can be done externally)
     * - in the learning mode, remove a target if too many are masked across the samples (in that case, max likelihood
     *   parameter estimation is unreliable)
     */
    private static ReadCountCollection processReadCountCollection(@Nonnull final ReadCountCollection rawReadCounts,
                                                                  @Nonnull final CoverageModelArgumentCollection args,
                                                                  @Nonnull final SexGenotypeDataCollection sexGenotypeDataCollection,
                                                                  @Nonnull final GermlinePloidyAnnotatedTargetCollection germlinePloidyAnnotatedTargetCollection) {
        /* stage 1. remove bad values */
        final ReadCountCollection readCountsStage_1 = ReadCountCollectionUtils.removeColumnsWithBadValues(rawReadCounts, logger);

        /* stage 2. remove targets that are masked for learning across all samples */
        final BiFunction<String, Target, Integer> germlinePloidyExtractorBiFunction = (sampleName, target) ->
                germlinePloidyAnnotatedTargetCollection.getTargetGermlinePloidyByGenotype(target,
                        sexGenotypeDataCollection.getSampleSexGenotypeData(sampleName).getSexGenotype());
        final BiFunction<Double, Integer, Boolean> isMaskedForLearningBiFunction = (readCount, germlinePloidy) ->
                isMaskedForLearning(readCount, germlinePloidy, args.getMinLearningReadCount());
        final ReadCountCollection readCountsStage_2 = removeMaskedForLearning(readCountsStage_1,
                germlinePloidyExtractorBiFunction, isMaskedForLearningBiFunction, logger);
        return readCountsStage_2;
    }

    /**
     * Removes targets that are masked for learning across all samples from a read count collection
     * @param readCounts original read counts
     * @param germlinePloidyExtractor (sampleName, Target) -> germline ploidy
     * @param isMaskedForLearning (read count, germline ploidy) -> boolean
     * @param logger a logger instance
     * @return filtered read count collection; may be a reference to the original collection
     */
    private static ReadCountCollection removeMaskedForLearning(@Nonnull final ReadCountCollection readCounts,
                                                               @Nonnull final BiFunction<String, Target, Integer> germlinePloidyExtractor,
                                                               @Nonnull final BiFunction<Double, Integer, Boolean> isMaskedForLearning,
                                                               @Nonnull final Logger logger) {
        final Predicate<ReadCountCollection.EntryIndex> isMaskedForLearningPredicate = entry -> isMaskedForLearning.apply(
                readCounts.counts().getEntry(entry.targetIndex, entry.sampleIndex),
                germlinePloidyExtractor.apply(readCounts.columnNames().get(entry.sampleIndex), readCounts.targets().get(entry.targetIndex)));
        return ReadCountCollectionUtils.removeTargetsWithTooManyFlags(readCounts, isMaskedForLearningPredicate,
                readCounts.counts().getColumnDimension() - 1, "masked for learning", logger);
    }

    /**
     * Partitions the target space into {@link #numTargetBlocks} contiguous blocks
     */
    private void initializeTargetSpaceBlocks() {
        targetBlocks = CoverageModelSparkUtils.createLinearlySpacedIndexBlocks(numTargets, numTargetBlocks,
                CoverageModelGlobalConstants.DEFAULT_MIN_TARGET_BLOCK_SIZE);
        logger.debug("Target space blocks: " + targetBlocks.stream().map(LinearlySpacedIndexBlock::toString)
                .reduce((L, R) -> L + "\t" + R).orElse("None"));
    }

    /**
     * Determines whether a target with read count {@code readCount} and germline ploidy {@code germlinePloidy}
     * should be masked in the parameter estimation ("learning") step
     *
     * @param readCount raw read count
     * @param germlinePloidy germline ploidy
     * @param minLearningReadCount minimum read count for considering the read count for learning
     * @return a boolean
     */
    private static boolean isMaskedForLearning(final double readCount, final int germlinePloidy,
                                               final int minLearningReadCount) {
        return readCount < minLearningReadCount || germlinePloidy == 0;
    }

    /**
     * Pushes the read count data, learning mask, mapping error, etc. to compute blocks
     *
     * @implNote This method is written in a seemingly unnatural way: read counts are first converted to a
     * list of {@link ReadCountRecord}, one record for each target that contains the read counts for
     * all included samples. The records are chopped into blocks in the target space and are pushed to
     * the compute blocks. This implementation makes it easy to bypass the read count collection in the
     * future and create the compute RDD directly from a combined read count table living on HDFS.
     */
    @UpdatesRDD
    private void pushInitialDataToComputeBlocks() {
        /* parse reads and initialize containers */
        logger.info("Collecting coverage data to the driver node...");
        final List<ReadCountRecord> recs = processedReadCounts.records();

        final List<Tuple2<LinearlySpacedIndexBlock, CoverageModelEMComputeBlock.InitialDataBlock>> dataBlockList =
                new ArrayList<>();

        targetBlockStream().forEach(tb -> {
            /* take a contiguous [targets in the block x all samples] chunk from the read count collection
             * and ravel it in Fortran order */
            final int[] rawReadCountBlock = IntStream.range(tb.getBegIndex(), tb.getEndIndex())
                    /* 1-to-STATE flat map of each target to the read counts of all samples */
                    .mapToObj(ti -> recs.get(ti).getDoubleCounts())
                    .flatMapToDouble(Arrays::stream)
                    .mapToInt(d -> (int)FastMath.round(d)) /* round to integer values */
                    .toArray();

            /* fetch the germline ploidy within the same contiguous block of reads */
            final int[] germlinePloidyBlock = IntStream.range(tb.getBegIndex(), tb.getEndIndex())
                    /* map target index to actual targets */
                    .mapToObj(processedTargetList::get)
                    /* 1-to-STATE flat map of each target to the germline ploidy of all samples */
                    .map(target -> processedSampleNameList.stream()
                            .map(sampleName -> sexGenotypeDataCollection
                                    .getSampleSexGenotypeData(sampleName).getSexGenotype())
                            .mapToInt(sampleSexGenotype -> germlinePloidyAnnotatedTargetCollection
                                    .getTargetGermlinePloidyByGenotype(target, sampleSexGenotype))
                            .toArray())
                    .flatMapToInt(Arrays::stream)
                    .toArray();

            /* If a target has zero germline ploidy, set its read counts to READ_COUNT_ON_ZERO_PLOIDY_TARGETS
             * regardless of what read count table provides (= mapping error) */
            IntStream.range(0, rawReadCountBlock.length)
                    .filter(idx -> germlinePloidyBlock[idx] == 0)
                    .forEach(idx -> rawReadCountBlock[idx] = CoverageModelGlobalConstants.READ_COUNT_ON_ZERO_PLOIDY_TARGETS);

            /* targets with low read count or zero ploidy are masked */
            final int[] maskBlock = IntStream.range(0, rawReadCountBlock.length)
                    .map(idx -> isMaskedForLearning(rawReadCountBlock[idx], germlinePloidyBlock[idx], config.getMinLearningReadCount()) ? 0 : 1)
                    .toArray();

            /* TODO github/gatk-protected issue #855 -- in the future, this must be replaced with a sample-
             * and target-specific value calculated from the mapping quality distribution of each target */
            final double[] mappingErrorRateBlock = IntStream.range(0, rawReadCountBlock.length)
                    .mapToDouble(idx -> config.getMappingErrorRate())
                    .toArray();

            /* add the block to list */
            dataBlockList.add(new Tuple2<>(tb, new CoverageModelEMComputeBlock.InitialDataBlock(rawReadCountBlock,
                    maskBlock, mappingErrorRateBlock)));
        });

        /* push to compute blocks */
        logger.info("Pushing coverage data to worker(s)...");
        joinWithWorkersAndMap(dataBlockList, p -> p._1.cloneWithInitializedData(p._2));
    }

    /**
     * Initialize the compute blocks with given model parameters
     *
     * Note:
     *
     * - If ARD is disabled, a null ARD coefficient INDArray will be passed to the workers
     * - If bias covariates are disabled, null mean and var bias covariates INDArray will be passed to the workers
     *
     * @param model coverage model parameters
     */
    private void initializeWorkersWithGivenModel(@Nonnull final CoverageModelParameters model) {
        logger.info("Pushing model parameters to worker(s)...");
        if (sparkContextIsAvailable) {
            initializeWorkersWithGivenModelSpark(model);
        } else {
            initializeWorkersWithGivenModelLocal(model);
        }
    }

    @UpdatesRDD
    private void initializeWorkersWithGivenModelSpark(@Nonnull final CoverageModelParameters model) {
        final Broadcast<CoverageModelParameters> broadcastedModel = ctx.broadcast(model);

        /* basic model nodes */
        mapWorkers(cb -> {
            final LinearlySpacedIndexBlock tb = cb.getTargetSpaceBlock();
            return cb
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t,
                            broadcastedModel.getValue().getTargetMeanBiasOnTargetBlock(tb))
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t,
                            broadcastedModel.getValue().getTargetUnexplainedVarianceOnTargetBlock(tb));
        });

        /* nodes relating to bias covariates */
        if (biasCovariatesEnabled) {
            mapWorkers(cb -> {
                final LinearlySpacedIndexBlock tb = cb.getTargetSpaceBlock();
                return cb
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl,
                                broadcastedModel.getValue().getMeanBiasCovariatesOnTargetBlock(tb));
            });
        }

        /* if ARD is enabled, inject ARD coefficients */
        if (ardEnabled) {
            mapWorkers(cb -> {
                final LinearlySpacedIndexBlock tb = cb.getTargetSpaceBlock();
                return cb
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.alpha_l,
                                broadcastedModel.getValue().getBiasCovariateARDCoefficients());
            });
        }
    }

    @UpdatesRDD
    private void initializeWorkersWithGivenModelLocal(@Nonnull final CoverageModelParameters model) {
        /* basic model nodes */
        mapWorkers(cb -> {
            final LinearlySpacedIndexBlock tb = cb.getTargetSpaceBlock();
            return cb
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t,
                            model.getTargetMeanBiasOnTargetBlock(tb))
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t,
                            model.getTargetUnexplainedVarianceOnTargetBlock(tb));
        });

        /* nodes relating to bias covariates */
        if (biasCovariatesEnabled) {
            mapWorkers(cb -> {
                final LinearlySpacedIndexBlock tb = cb.getTargetSpaceBlock();
                return cb
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl,
                                model.getMeanBiasCovariatesOnTargetBlock(tb));
            });
        }

        /* if ARD is enabled, inject ARD coefficients and the covariance of bias covariates from the model as well */
        if (ardEnabled) {
            mapWorkers(cb -> cb
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.alpha_l,
                            model.getBiasCovariateARDCoefficients()));
        }
    }

    @UpdatesRDD
    private void initializeWorkerModelParameters() {
        if (processedModel != null) { /* an initial model is already provided */
            logger.info("Initializing workers with the provided model...");
            initializeWorkersWithGivenModel(processedModel);
        } else { /* model parameters need to be learned */
            Utils.validateArg(numLatents < numSamples, "Number of bias latent variables must be strictly less than the" +
                    " number of samples");
            if (config.getModelInitializationStrategy().equals(CoverageModelArgumentCollection.ModelInitializationStrategy.RANDOM) ||
                    !biasCovariatesEnabled) {
                initializeWorkersWithRandomModel();
            } else if (config.getModelInitializationStrategy().equals(CoverageModelArgumentCollection.ModelInitializationStrategy.PCA)) {
                initializeWorkersWithPCA();
            } else {
                throw new GATKException.ShouldNeverReachHereException("Bad value of model initialization strategy");
            }
        }
    }

    @UpdatesRDD
    private void initializeWorkerPosteriors() {
        /* make a local copy for lambda capture (these are small, so no broadcasting is necessary) */
        final INDArray sampleMeanLogReadDepths = this.sampleMeanLogReadDepths;
        final INDArray sampleVarLogReadDepths = this.sampleVarLogReadDepths;
        final INDArray sampleUnexplainedVariance = this.sampleUnexplainedVariance;
        final INDArray sampleBiasLatentPosteriorFirstMoments = this.sampleBiasLatentPosteriorFirstMoments;
        final INDArray sampleBiasLatentPosteriorSecondMoments = this.sampleBiasLatentPosteriorSecondMoments;

        /* calculate copy ratio prior expectations */
        logger.info("Calculating copy ratio priors on the driver node...");
        final List<CopyRatioExpectations> copyRatioPriorExpectationsList = sampleIndexStream()
                .mapToObj(si -> copyRatioExpectationsCalculator.getCopyRatioPriorExpectations(
                        CopyRatioCallingMetadata.builder()
                                .sampleName(processedSampleNameList.get(si))
                                .sampleSexGenotypeData(processedSampleSexGenotypeData.get(si))
                                .sampleAverageMappingErrorProbability(config.getMappingErrorRate())
                                .build(), processedTargetList))
                .collect(Collectors.toList());

        /* update log chain posterior expectation */
        sampleLogChainPosteriors.assign(Nd4j.create(copyRatioPriorExpectationsList.stream()
                .mapToDouble(CopyRatioExpectations::getLogChainPosteriorProbability)
                .toArray(), new int[] {numSamples, 1}));

        /* push per-target copy ratio expectations to workers */
        final List<Tuple2<LinearlySpacedIndexBlock, Tuple2<INDArray, INDArray>>> copyRatioPriorsList = new ArrayList<>();
        for (final LinearlySpacedIndexBlock tb : targetBlocks) {
            final double[] logCopyRatioPriorMeansBlock = IntStream.range(0, tb.getNumElements())
                    .mapToObj(rti -> copyRatioPriorExpectationsList.stream()
                            .mapToDouble(cre -> cre.getLogCopyRatioMeans()[rti + tb.getBegIndex()])
                            .toArray())
                    .flatMapToDouble(Arrays::stream)
                    .toArray();

            final double[] logCopyRatioPriorVariancesBlock = IntStream.range(0, tb.getNumElements())
                    .mapToObj(rti -> copyRatioPriorExpectationsList.stream()
                            .mapToDouble(cre -> cre.getLogCopyRatioVariances()[rti + tb.getBegIndex()])
                            .toArray())
                    .flatMapToDouble(Arrays::stream)
                    .toArray();

            /* we do not need to take care of log copy ratio means and variances on masked targets here.
             * potential NaNs will be rectified in the compute blocks by calling the method
             * {@link CoverageModelEMComputeBlock#cloneWithInitializedData} */

            copyRatioPriorsList.add(new Tuple2<>(tb, new Tuple2<INDArray, INDArray>(
                    Nd4j.create(logCopyRatioPriorMeansBlock, new int[] {numSamples, tb.getNumElements()}, 'f'),
                    Nd4j.create(logCopyRatioPriorVariancesBlock, new int[] {numSamples, tb.getNumElements()}, 'f'))));
        }

        /* push to compute blocks */
        logger.info("Pushing posteriors to worker(s)...");

        /* copy ratio priors */
        joinWithWorkersAndMap(copyRatioPriorsList, p -> p._1.cloneWithUpdateCopyRatioPriors(p._2._1, p._2._2));

        /* read depth and sample-specific unexplained variance */
        mapWorkers(cb -> cb
                .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.log_d_s,
                        sampleMeanLogReadDepths)
                .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.var_log_d_s,
                        sampleVarLogReadDepths)
                .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.gamma_s,
                        sampleUnexplainedVariance));

        /* if bias covariates are enabled, initialize E[z] and E[z z^T] as well */
        if (biasCovariatesEnabled) {
            mapWorkers(cb -> cb
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.z_sl,
                            sampleBiasLatentPosteriorFirstMoments)
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.zz_sll,
                            sampleBiasLatentPosteriorSecondMoments));
        }
    }

    /**
     * Initialize workers with randomly generated model parameters.
     */
    private void initializeWorkersWithRandomModel() {
        logger.info("Initializing workers with a random model...");
        final CoverageModelParameters randomModel = CoverageModelParameters.generateRandomModel(
                processedTargetList, numLatents,
                CoverageModelGlobalConstants.RANDOM_MODEL_SEED,
                CoverageModelGlobalConstants.RANDOM_MEAN_LOG_BIAS_STD,
                CoverageModelGlobalConstants.RANDOM_BIAS_COVARIATES_STD,
                CoverageModelGlobalConstants.RANDOM_UNEXPLAINED_VARIANCE_MAX,
                biasCovariatesARDCoefficients);
        initializeWorkersWithGivenModel(randomModel);
    }

    /**
     * Initialize model parameters by performing PCA.
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    private void initializeWorkersWithPCA() {
        logger.info("Initializing model parameters using PCA...");
        /* initially, set m_t, Psi_t and W_tl to zero to get an estimate of the read depth */
        final int numLatents = config.getNumLatents();
        mapWorkers(cb -> cb
                .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t,
                        Nd4j.zeros(new int[] {1, cb.getTargetSpaceBlock().getNumElements()}))
                .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t,
                        Nd4j.zeros(new int[] {1, cb.getTargetSpaceBlock().getNumElements()})));
        if (biasCovariatesEnabled) {
            mapWorkers(cb -> cb
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl,
                            Nd4j.zeros(new int[] {cb.getTargetSpaceBlock().getNumElements(), numLatents})));
        }

        /* update read depth without taking into account correction from bias covariates */
        updateReadDepthPosteriorExpectations(1.0, true);

        /* fetch sample covariance matrix */
        final int minPCAInitializationReadCount = config.getMinPCAInitializationReadCount();
        mapWorkers(cb -> cb.cloneWithPCAInitializationData(minPCAInitializationReadCount, Integer.MAX_VALUE));
        cacheWorkers("PCA initialization");
        final INDArray targetCovarianceMatrix = mapWorkersAndReduce(
                CoverageModelEMComputeBlock::calculateTargetCovarianceMatrixForPCAInitialization,
                INDArray::add);

        /* perform eigen-decomposition on the target covariance matrix */
        final ImmutablePair<INDArray, INDArray> targetCovarianceEigensystem = CoverageModelEMWorkspaceMathUtils.eig(
                targetCovarianceMatrix, false, logger);

        /* the eigenvalues of sample covariance matrix can be immediately inferred by scaling */
        final INDArray sampleCovarianceEigenvalues = targetCovarianceEigensystem.getLeft().div(numSamples);

        /* estimate the isotropic unexplained variance -- see Bishop 12.46 */
        final int residualDim = numTargets - numLatents;
        final double isotropicVariance = sampleCovarianceEigenvalues.get(NDArrayIndex.interval(numLatents, numSamples))
                .sumNumber().doubleValue() / residualDim;
        logger.info(String.format("PCA estimate of isotropic unexplained variance: %f" , isotropicVariance));

        /* estimate bias factors -- see Bishop 12.45 */
        final INDArray scaleFactors = Transforms.sqrt(sampleCovarianceEigenvalues
                .get(NDArrayIndex.interval(0, numLatents)).sub(isotropicVariance), false);
        final INDArray biasCovariatesPCA = Nd4j.create(new int[] {numTargets, numLatents});
        for (int li = 0; li < numLatents; li++) {
            final INDArray v = targetCovarianceEigensystem.getRight().getColumn(li);
            /* calculate [Delta_PCA_st]^T v */
            /* note: we do not need to broadcast vec since it is small and lambda capture is just fine */
            final INDArray unnormedBiasCovariate = CoverageModelSparkUtils.assembleINDArrayBlocksFromCollection(
                    mapWorkersAndCollect(cb -> ImmutablePair.of(cb.getTargetSpaceBlock(),
                            cb.getINDArrayFromCache(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Delta_PCA_st)
                                    .transpose().mmul(v))), 0);
            final double norm = unnormedBiasCovariate.norm1Number().doubleValue();
            final INDArray normedBiasCovariate = unnormedBiasCovariate
                    .divi(norm)
                    .muli(scaleFactors.getDouble(li));
            biasCovariatesPCA.getColumn(li).assign(normedBiasCovariate);
        }
        if (ardEnabled) { /* a better estimate of ARD coefficients */
            biasCovariatesARDCoefficients.assign(Nd4j.zeros(new int[]{1, numLatents})
                    .addi(config.getInitialARDPrecisionRelativeToNoise() / isotropicVariance));
        }

        final CoverageModelParameters modelParamsFromPCA = new CoverageModelParameters(
                processedTargetList,
                Nd4j.zeros(new int[] {1, numTargets}),
                Nd4j.zeros(new int[] {1, numTargets}).addi(isotropicVariance),
                biasCovariatesPCA,
                biasCovariatesARDCoefficients);

        /* clear PCA initialization data from workers */
        mapWorkers(CoverageModelEMComputeBlock::cloneWithRemovedPCAInitializationData);

        /* push model parameters to workers */
        initializeWorkersWithGivenModel(modelParamsFromPCA);

        /* update bias latent posterior expectations without admixing */
        updateBiasLatentPosteriorExpectations(1.0);
    }

    /**
     * This method fetches bias covariates from the compute blocks, applies the Fourier filter on it,
     * partitions the result in the target space, and pushes it to compute block(s)
     */
    @EvaluatesRDD @UpdatesRDD
    private void updateFilteredBiasCovariates() {
        updateFilteredBiasCovariates(fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl, 0));
    }

    /**
     * This method applies the Fourier filter on a given bias covariates matrix, applies the Fourier filter on it,
     * partitions the result, and pushes it to compute block(s)
     *
     * @param biasCovariates any T x D bias covariates matrix
     */
    @UpdatesRDD
    private void updateFilteredBiasCovariates(@Nonnull final INDArray biasCovariates) {
        final INDArray filteredBiasCovariates = Nd4j.create(biasCovariates.shape());

        /* instantiate the Fourier filter */
        final FourierLinearOperatorNDArray regularizerFourierLinearOperator = createRegularizerFourierLinearOperator();

        /* FFT by resolving W_tl on l */
        for (int li = 0; li < numLatents; li++) {
            final INDArrayIndex[] slice = {NDArrayIndex.all(), NDArrayIndex.point(li)};
            filteredBiasCovariates.get(slice).assign(
                    regularizerFourierLinearOperator.operate(biasCovariates.get(slice)));
        }

        /* sent the new W to workers */
        switch (config.getBiasCovariatesComputeNodeCommunicationPolicy()) {
            case BROADCAST_HASH_JOIN:
                pushToWorkers(mapINDArrayToBlocks(filteredBiasCovariates),
                        (W, cb) -> cb.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock
                                .CoverageModelICGCacheNode.F_W_tl, W.get(cb.getTargetSpaceBlock())));
                break;

            case RDD_JOIN:
                joinWithWorkersAndMap(chopINDArrayToBlocks(filteredBiasCovariates),
                        p -> p._1.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock
                                .CoverageModelICGCacheNode.F_W_tl, p._2));
                break;
        }
    }

    /**
     * Creates the Fourier operator for CNV-avoiding regularizier
     * @return an instance of {@link FourierLinearOperatorNDArray}
     */
    @EvaluatesRDD
    private FourierLinearOperatorNDArray createRegularizerFourierLinearOperator() {
        final double psiAverage = mapWorkersAndReduce(cb -> cb.getINDArrayFromCache(CoverageModelEMComputeBlock
                .CoverageModelICGCacheNode.M_Psi_inv_st)
                .sumNumber().doubleValue(), (d1, d2) -> d1 + d2) / (numTargets * numSamples);
        final double fact = config.getFourierRegularizationStrength() * psiAverage;
        return new FourierLinearOperatorNDArray(numTargets,
                Arrays.stream(fourierFactors).map(f -> f * fact).toArray(), config.zeroPadFFT());
    }

    /* E-step methods */

    /**
     * Updates the first (E[z]) and second (E[z z^T]) posterior moments of the log bias continuous
     * latent variables (z)
     *
     * @implNote the operations done on the driver node have low complexity only if D, the dimension of the latent
     * space, is small:
     *
     *     (a) G_s = (I + [G_partial_sll])^{-1} for each sample \sim O(STATE x D^3)
     *     (b) E[z_s] = G_s [z_rhs_ls] for each sample \sim O(STATE x D^3)
     *     (c) E[z_s z_s^T] = G_s + E[z_s] E[z_s^T] for each sample \sim O(STATE x D^2)
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateBiasLatentPosteriorExpectations(final double admixingRatio) {
        /* calculate and cache the required quantities on compute blocks */
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_Z));
        cacheWorkers("after E-step for bias initialization");

        final ImmutableTriple<INDArray, INDArray, INDArray> biasLatentPosteriorExpectationsData =
                fetchBiasLatentPosteriorExpectationsDataFromWorkers();
        final INDArray G_partial_sll = biasLatentPosteriorExpectationsData.left;
        final INDArray z_rhs_ls = biasLatentPosteriorExpectationsData.middle;
        final INDArray shared_ll = biasLatentPosteriorExpectationsData.right;

        /* calculate G_{s\mu\nu} = (sharedPart + W^T M \Psi^{-1} W)^{-1} by doing sample-wise matrix inversion */
        final INDArray new_z_sl = Nd4j.create(new int[] {numSamples, numLatents});
        final INDArray new_zz_sll = Nd4j.create(new int[] {numSamples, numLatents, numLatents});

        for (int si = 0; si < numSamples; si++) {
            final INDArray G_ll = CoverageModelEMWorkspaceMathUtils.minv(shared_ll
                    .add(G_partial_sll.get(NDArrayIndex.point(si), NDArrayIndex.all(), NDArrayIndex.all())));
            /* E[z_s] = G_s W^T M_{st} \Psi_{st}^{-1} (m_{st} - m_t) */
            new_z_sl.getRow(si).assign(G_ll.mmul(z_rhs_ls.getColumn(si)).transpose());
            /* E[z_s z_s^T] = G_s + E[z_s] E[z_s^T] */
            final INDArray z = new_z_sl.get(NDArrayIndex.point(si), NDArrayIndex.all());
            new_zz_sll.get(NDArrayIndex.point(si), NDArrayIndex.all(), NDArrayIndex.all())
                    .assign(G_ll.add(z.transpose().mmul(z)));
        }

        /* admix with old posteriors */
        final INDArray new_z_sl_admixed = new_z_sl
                .mul(admixingRatio).addi(sampleBiasLatentPosteriorFirstMoments.mul(1.0 - admixingRatio));
        final INDArray new_zz_sll_admixed = new_zz_sll
                .mul(admixingRatio).addi(sampleBiasLatentPosteriorSecondMoments.mul(1.0 - admixingRatio));

        /* calculate the error from the change in E[z_s] */
        final double errorNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                new_z_sl_admixed.sub(sampleBiasLatentPosteriorFirstMoments));

        /* update driver-node copies */
        sampleBiasLatentPosteriorFirstMoments.assign(new_z_sl_admixed);
        sampleBiasLatentPosteriorSecondMoments.assign(new_zz_sll_admixed);

        /* broadcast the new bias latent posteriors */
        pushToWorkers(ImmutablePair.of(new_z_sl_admixed, new_zz_sll_admixed),
                (dat, cb) -> cb
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.z_sl, dat.left)
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.zz_sll, dat.right));

        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .build();
    }

    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateBiasLatentPosteriorExpectations() {
        return updateBiasLatentPosteriorExpectations(config.getMeanFieldAdmixingRatio());
    }

    /**
     * Fetches the data required for calculating bias latent posterior expectations from the compute blocks
     * The return value is an {@link ImmutableTriple} of (contribGMatrix, contribZ, sharedPart)
     *
     * For the definition of the first two elements, refer to the javadoc of
     * {@link CoverageModelEMComputeBlock#getBiasLatentPosteriorDataUnregularized}
     *
     * If the Fourier regularizer is enabled, sharedPart = [I] + W^T \lambda [F] [W]
     * If the Fourier regularizer is disabled, sharedPart = [I]
     *
     * @return an {@link ImmutableTriple}
     */
    private ImmutableTriple<INDArray, INDArray, INDArray> fetchBiasLatentPosteriorExpectationsDataFromWorkers() {
        /* query the compute blocks for bias latent posterior data and reduce by pairwise addition */
        if (config.fourierRegularizationEnabled()) {
            /* TODO github/gatk-protected issue #853 -- it is desirable to perform in-place reduction but Spark
             * looses resilience; perhaps we can use accumulators instead of map-reduce */
            final ImmutableTriple<INDArray, INDArray, INDArray> data =
                    mapWorkersAndReduce(CoverageModelEMComputeBlock::getBiasLatentPosteriorDataRegularized,
                            (t1, t2) -> ImmutableTriple.of(t1.left.add(t2.left), t1.middle.add(t2.middle),
                                    t1.right.add(t2.right)));
            /* sharedPart = I + W^T [\lambda F W] */
            final INDArray contribFilter = data.right;
            final INDArray sharedPart = Nd4j.eye(numLatents).addi(contribFilter);
            return ImmutableTriple.of(data.left, data.middle, sharedPart);
        } else { /* w/o regularizer */
            final ImmutablePair<INDArray, INDArray> data =
                    mapWorkersAndReduce(CoverageModelEMComputeBlock::getBiasLatentPosteriorDataUnregularized,
                            (p1, p2) -> ImmutablePair.of(p1.left.add(p2.left), p1.right.add(p2.right)));
            /* sharedPart = I */
            final INDArray sharedPart = Nd4j.eye(numLatents);
            return ImmutableTriple.of(data.left, data.right, sharedPart);
        }
    }

    /**
     * E-step update of read depth ($d_s$)
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateReadDepthPosteriorExpectations(final double admixingRatio,
                                                                 final boolean neglectBiasCovariates) {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_D));
        cacheWorkers("after E-step for read depth initialization");
        /* map each compute block to their respective read depth estimation data (a triple of rank-1 INDArray's each
         * with STATE elements) and reduce by pairwise addition */
        final ImmutablePair<INDArray, INDArray> factors = mapWorkersAndReduce(
                cb -> cb.getReadDepthLatentPosteriorData(neglectBiasCovariates),
                (p1, p2) -> ImmutablePair.of(p1.left.add(p2.left), p1.right.add(p2.right)));

        /* put together */
        final INDArray numerator = factors.left;
        final INDArray denominator = factors.right;

        final INDArray newSampleMeanLogReadDepths = numerator.div(denominator);
        final INDArray newSampleVarLogReadDepths = Nd4j.ones(denominator.shape()).div(denominator);

        final INDArray newSampleMeanLogReadDepthsAdmixed;
        final INDArray newSampleVarLogReadDepthsAdmixed;

        /* admix */
        newSampleMeanLogReadDepthsAdmixed = newSampleMeanLogReadDepths
                .mul(admixingRatio)
                .addi(sampleMeanLogReadDepths.mul(1.0 - admixingRatio));
        newSampleVarLogReadDepthsAdmixed = newSampleVarLogReadDepths
                .mul(admixingRatio)
                .addi(sampleVarLogReadDepths.mul(1.0 - admixingRatio));

        /* calculate the error only using the change in the mean */
        final double errorNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                newSampleMeanLogReadDepthsAdmixed.sub(sampleMeanLogReadDepths));

        /* update local copies of E[log(d_s)] and var[log(d_s)] */
        sampleMeanLogReadDepths.assign(newSampleMeanLogReadDepthsAdmixed);
        sampleVarLogReadDepths.assign(newSampleVarLogReadDepthsAdmixed);

        /* push E[log(d_s)] and var[log(d_s)] to all workers; they all need a copy */
        pushToWorkers(ImmutablePair.of(newSampleMeanLogReadDepths, newSampleVarLogReadDepths),
                (dat, cb) -> cb
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.log_d_s, dat.left)
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.var_log_d_s, dat.right));

        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .build();
    }

    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateReadDepthPosteriorExpectations() {
        return updateReadDepthPosteriorExpectations(config.getMeanFieldAdmixingRatio(), false);
    }

    /**
     * E-step update of the sample-specific unexplained variance
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm") and the average
     * number of function evaluations per sample (key: "iterations")
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateSampleUnexplainedVariance() {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_GAMMA));
        cacheWorkers("after E-step for sample unexplained variance initialization");

        /* create a compound objective function for simultaneous multi-sample queries */
        final java.util.function.Function<Map<Integer, Double>, Map<Integer, Double>> objFunc = arg -> {
            if (arg.isEmpty()) { /* nothing to evaluate */
                return Collections.emptyMap();
            }
            /* get sample indices */
            final int[] sampleIndices = arg.keySet().stream().mapToInt(i -> i).toArray();
            /* get values of gamma for each sample index */
            final INDArray gammaValues = Nd4j.create(Arrays.stream(sampleIndices)
                    .mapToDouble(arg::get).toArray(), new int[] {sampleIndices.length, 1});
            /* query */
            final INDArray eval = mapWorkersAndReduce(cb ->
                    cb.calculateSampleSpecificVarianceObjectiveFunctionMultiSample(sampleIndices, gammaValues), INDArray::add);
            /* create output map */
            final Map<Integer, Double> output = new HashMap<>();
            IntStream.range(0, sampleIndices.length)
                    .forEach(evalIdx -> output.put(sampleIndices[evalIdx], eval.getDouble(evalIdx)));
            return output;
        };

        final java.util.function.Function<UnivariateSolverSpecifications, AbstractUnivariateSolver> solverFactory = spec ->
                new RobustBrentSolver(spec.getRelativeAccuracy(), spec.getAbsoluteAccuracy(),
                        spec.getFunctionValueAccuracy(), null, config.getSampleSpecificVarianceSolverNumBisections(),
                        config.getSampleSpecificVarianceSolverRefinementDepth());

        /* instantiate a synchronized multi-sample root finder and add jobs */
        final SynchronizedUnivariateSolver syncSolver = new SynchronizedUnivariateSolver(objFunc, solverFactory, numSamples);
        IntStream.range(0, numSamples)
                .forEach(si -> {
                    final double x0 = 0.5 * config.getSampleSpecificVarianceUpperLimit();
                    syncSolver.add(si, 0, config.getSampleSpecificVarianceUpperLimit(), x0,
                            config.getSampleSpecificVarianceAbsoluteTolerance(), config.getSampleSpecificVarianceRelativeTolerance(),
                            config.getSampleSpecificVarianceMaximumIterations());
                });

        /* solve and collect statistics */
        final INDArray newSampleUnexplainedVariance = Nd4j.create(numSamples, 1);
        final List<Integer> numberOfEvaluations = new ArrayList<>(numSamples);
        try {
            final Map<Integer, SynchronizedUnivariateSolver.UnivariateSolverSummary> newSampleSpecificVarianceMap = syncSolver.solve();
            newSampleSpecificVarianceMap.entrySet().forEach(entry -> {
                final int sampleIndex = entry.getKey();
                final SynchronizedUnivariateSolver.UnivariateSolverSummary summary = entry.getValue();
                double val =  0;
                switch (summary.status) {
                    case SUCCESS:
                        val = summary.x;
                        break;
                    case TOO_MANY_EVALUATIONS:
                        logger.warn("Could not locate the root of gamma -- increase the maximum number of" +
                                "function evaluations");
                        break;
                }
                newSampleUnexplainedVariance.put(sampleIndex, 0, val);
                numberOfEvaluations.add(summary.evaluations);
            });
        } catch (final InterruptedException ex) {
            throw new RuntimeException("The update of sample unexplained variance was interrupted -- can not continue");
        }

        /* admix */
        final INDArray newSampleUnexplainedVarianceAdmixed = newSampleUnexplainedVariance
                .mul(config.getMeanFieldAdmixingRatio())
                .addi(sampleUnexplainedVariance.mul(1 - config.getMeanFieldAdmixingRatio()));

        /* calculate the error */
        final double errorNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                newSampleUnexplainedVarianceAdmixed.sub(sampleUnexplainedVariance));

        /* update local copy */
        sampleUnexplainedVariance.assign(newSampleUnexplainedVarianceAdmixed);

        /* push to workers */
        pushToWorkers(newSampleUnexplainedVarianceAdmixed, (arr, cb) -> cb.cloneWithUpdatedPrimitive(
                CoverageModelEMComputeBlock.CoverageModelICGCacheNode.gamma_s,
                newSampleUnexplainedVarianceAdmixed));

        final int iterations = (int)(numberOfEvaluations.stream().mapToDouble(d -> d).sum() / numSamples);
        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .put(StandardSubroutineSignals.ITERATIONS, iterations)
                .build();
    }

    /**
     * E-step update of copy ratio posteriors
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateCopyRatioPosteriorExpectations(final double admixingRatio) {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_C));
        cacheWorkers("after E-step for copy ratio initialization");

        /* calculate posteriors */
        long startTime = System.nanoTime();
        final SubroutineSignal sig;
        if (config.getCopyRatioHMMType().equals(CoverageModelArgumentCollection.CopyRatioHMMType.LOCAL) || !sparkContextIsAvailable) {
            /* local mode */
            sig = updateCopyRatioPosteriorExpectationsLocal(admixingRatio);
        } else {
            /* spark mode */
            sig = updateCopyRatioPosteriorExpectationsSpark(admixingRatio);
        }
        long endTime = System.nanoTime();
        logger.debug("Copy ratio posteriors calculation time: " + (double)(endTime - startTime)/1000000 + " ms");
        return sig;
    }

    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateCopyRatioPosteriorExpectations() {
        return updateCopyRatioPosteriorExpectations(config.getMeanFieldAdmixingRatio());
    }

    /**
     * Fetches forward-backward and Viterbi algorithm results on all samples
     *
     * @return a list of {@link CopyRatioHMMResults}
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    protected List<List<HiddenStateSegmentRecord<STATE, Target>>> getCopyRatioSegments() {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_C));
        cacheWorkers("after E-step for copy ratio segment generation");
        final List<List<HiddenStateSegmentRecord<STATE, Target>>> result;
        /* calculate posteriors */
        long startTime = System.nanoTime();
        if (config.getCopyRatioHMMType().equals(CoverageModelArgumentCollection.CopyRatioHMMType.LOCAL) || !sparkContextIsAvailable) {
            /* local mode */
            result = getCopyRatioSegmentsLocal();
        } else {
            /* spark mode */
            result = getCopyRatioSegmentsSpark();
        }
        long endTime = System.nanoTime();
        logger.debug("Copy ratio segmentation time: " + (double)(endTime - startTime)/1000000 + " ms");
        return result;
    }

    /**
     * Local implementation of the E-step update of copy ratio posteriors
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    public SubroutineSignal updateCopyRatioPosteriorExpectationsLocal(final double admixingRatio) {
        /* step 1. fetch copy ratio emission data */
        final List<List<CoverageModelCopyRatioEmissionData>> copyRatioEmissionData = fetchCopyRatioEmissionDataLocal();

        /* step 2. run the forward-backward algorithm and calculate copy ratio posteriors */
        final INDArray sampleReadDepths = Transforms.exp(sampleMeanLogReadDepths, true);
        final List<CopyRatioExpectations> copyRatioPosteriorResults = sampleIndexStream()
                .parallel()
                .mapToObj(si -> copyRatioExpectationsCalculator.getCopyRatioPosteriorExpectations(
                        CopyRatioCallingMetadata.builder()
                                .sampleName(processedSampleNameList.get(si))
                                .sampleSexGenotypeData(processedSampleSexGenotypeData.get(si))
                                .sampleCoverageDepth(sampleReadDepths.getDouble(si))
                                .emissionCalculationStrategy(EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                                .build(),
                        processedTargetList,
                        copyRatioEmissionData.get(si)))
                .collect(Collectors.toList());

        /* update log chain posterior expectation */
        sampleLogChainPosteriors.assign(Nd4j.create(copyRatioPosteriorResults.stream()
                .mapToDouble(CopyRatioExpectations::getLogChainPosteriorProbability)
                .toArray(), new int[] {numSamples, 1}));

        /* sent the results back to workers */
        final ImmutablePair<INDArray, INDArray> copyRatioPosteriorDataPair =
                convertCopyRatioLatentPosteriorExpectationsToNDArray(copyRatioPosteriorResults);
        final INDArray log_c_st = copyRatioPosteriorDataPair.left;
        final INDArray var_log_c_st = copyRatioPosteriorDataPair.right;

        /* partition the pair of (log_c_st, var_log_c_st), sent the result to workers via broadcast-hash-map */
        pushToWorkers(mapINDArrayPairToBlocks(log_c_st.transpose(), var_log_c_st.transpose()),
                (p, cb) -> cb.cloneWithUpdatedCopyRatioPosteriors(
                        p.get(cb.getTargetSpaceBlock()).left.transpose(),
                        p.get(cb.getTargetSpaceBlock()).right.transpose(),
                        admixingRatio));
        cacheWorkers("after E-step update of copy ratio posteriors");

        /* collect subroutine signals */
        final List<SubroutineSignal> sigs = mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal);

        final double errorNormInfinity = Collections.max(sigs.stream()
                .map(sig -> sig.<Double>get(StandardSubroutineSignals.RESIDUAL_ERROR_NORM))
                .collect(Collectors.toList()));

        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .build();
    }

    /**
     * Queries copy ratio emission data from compute blocks
     *
     * @return a double list of {@link CoverageModelCopyRatioEmissionData}
     */
    private List<List<CoverageModelCopyRatioEmissionData>> fetchCopyRatioEmissionDataLocal() {
        /* fetch data from workers */
        final List<ImmutablePair<LinearlySpacedIndexBlock, List<List<CoverageModelCopyRatioEmissionData>>>> collectedCopyRatioData =
                mapWorkersAndCollect(cb -> ImmutablePair.of(cb.getTargetSpaceBlock(), cb.getSampleCopyRatioLatentPosteriorData()));
        /* assemble and return */
        return IntStream.range(0, numSamples).parallel()
                .mapToObj(si -> CoverageModelSparkUtils.assembleBlockifiedCollection(
                        collectedCopyRatioData.stream()
                                .map(p -> ImmutablePair.of(p.getKey(), p.getValue().get(si)))
                                .collect(Collectors.toList())))
                .collect(Collectors.toList());
    }

    /**
     * Converts a list of copy ratio posteriors into a pair of INDArrays (log_c_st, var_log_c_st)
     *
     * @param copyRatioExpectationsList a list of {@link CopyRatioExpectations}
     * @return a pair of INDArrays (log_c_st, var_log_c_st)
     */
    private ImmutablePair<INDArray, INDArray> convertCopyRatioLatentPosteriorExpectationsToNDArray(
            @Nonnull final List<CopyRatioExpectations> copyRatioExpectationsList) {

        final INDArray log_c_st = Nd4j.create(numSamples, numTargets);
        final INDArray var_log_c_st = Nd4j.create(numSamples, numTargets);

        for (int si = 0; si < numSamples; si++) {
            final CopyRatioExpectations res = copyRatioExpectationsList.get(si);
            log_c_st.getRow(si).assign(Nd4j.create(res.getLogCopyRatioMeans(), new int[] {1, numTargets}));
            var_log_c_st.getRow(si).assign(Nd4j.create(res.getLogCopyRatioVariances(), new int[] {1, numTargets}));
        }
        return ImmutablePair.of(log_c_st, var_log_c_st);
    }

    private List<List<HiddenStateSegmentRecord<STATE, Target>>> getCopyRatioSegmentsLocal() {
        final List<List<CoverageModelCopyRatioEmissionData>> copyRatioEmissionData = fetchCopyRatioEmissionDataLocal();
        final INDArray sampleReadDepths = Transforms.exp(sampleMeanLogReadDepths, true);
        return sampleIndexStream()
                .mapToObj(si -> {
                    final CopyRatioCallingMetadata metadata = CopyRatioCallingMetadata.builder()
                            .sampleName(processedSampleNameList.get(si))
                            .sampleSexGenotypeData(processedSampleSexGenotypeData.get(si))
                            .sampleCoverageDepth(sampleReadDepths.getDouble(si))
                            .emissionCalculationStrategy(EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                            .build();
                    return copyRatioExpectationsCalculator.getCopyRatioHMMResults(metadata,
                            processedTargetList, copyRatioEmissionData.get(si));
                })
                /* segment each sample individually */
                .map(result -> {
                    final HMMSegmentProcessor<CoverageModelCopyRatioEmissionData, STATE, Target> processor =
                            new HMMSegmentProcessor<>(
                                    Collections.singletonList(result.getMetaData().getSampleName()),
                                    Collections.singletonList(result.getMetaData().getSampleSexGenotypeData()),
                                    referenceStateFactory,
                                    Collections.singletonList(new HashedListTargetCollection<>(processedTargetList)),
                                    Collections.singletonList(result.getForwardBackwardResult()),
                                    Collections.singletonList(result.getViterbiResult()));
                    return processor.getSegmentsAsList();
                })
                .collect(Collectors.toList());
    }

    /**
     * The Spark implementation of the E-step update of copy ratio posteriors
     *
     * @return a {@link SubroutineSignal} containing the update size
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    private SubroutineSignal updateCopyRatioPosteriorExpectationsSpark(final double admixingRatio) {
        /* local final member variables for lambda capture */
        final List<LinearlySpacedIndexBlock> targetBlocks = new ArrayList<>();
        targetBlocks.addAll(this.targetBlocks);
        final List<Target> targetList = new ArrayList<>();
        targetList.addAll(processedTargetList);
        final List<String> sampleNameList = new ArrayList<>();
        sampleNameList.addAll(processedSampleNameList);
        final List<SexGenotypeData> sampleSexGenotypeData = new ArrayList<>();
        sampleSexGenotypeData.addAll(processedSampleSexGenotypeData);
        final int numTargetBlocks = targetBlocks.size();
        final CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, STATE> calculator =
                this.copyRatioExpectationsCalculator;
        final INDArray sampleReadDepths = Transforms.exp(sampleMeanLogReadDepths, true);

        /* make an RDD of copy ratio posterior expectations */
        final JavaPairRDD<Integer, CopyRatioExpectations> copyRatioPosteriorExpectationsPairRDD =
                /* fetch copy ratio emission data from workers */
                fetchCopyRatioEmissionDataSpark()
                        /* calculate copy ratio posterior expectations; the original partitioning is preserved
                         * in the mapPartitions operation */
                        .mapPartitionsToPair(it -> {
                            final List<Tuple2<Integer, CopyRatioExpectations>> newPartitionData = new ArrayList<>();
                            while (it.hasNext()) {
                                final Tuple2<Integer, List<CoverageModelCopyRatioEmissionData>> prevDatum = it.next();
                                final int si = prevDatum._1;
                                final CopyRatioCallingMetadata copyRatioCallingMetadata = CopyRatioCallingMetadata.builder()
                                        .sampleName(sampleNameList.get(si))
                                        .sampleSexGenotypeData(sampleSexGenotypeData.get(si))
                                        .sampleCoverageDepth(sampleReadDepths.getDouble(si))
                                        .emissionCalculationStrategy(EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                                        .build();
                                newPartitionData.add(new Tuple2<>(prevDatum._1,
                                        calculator.getCopyRatioPosteriorExpectations(copyRatioCallingMetadata,
                                                targetList, prevDatum._2)));
                            }
                            return newPartitionData.iterator();
                        }, true);

        /* we need to do two things to copyRatioPosteriorExpectationsPairRDD; so we cache it */
        /* step 1. update log chain posterior expectation on the driver node */
        final double[] newSampleLogChainPosteriors = copyRatioPosteriorExpectationsPairRDD
                .mapValues(CopyRatioExpectations::getLogChainPosteriorProbability)
                .collect()
                .stream()
                .sorted(Comparator.comparingInt(t -> t._1))
                .mapToDouble(t -> t._2)
                .toArray();
        sampleLogChainPosteriors.assign(Nd4j.create(newSampleLogChainPosteriors, new int[] {numSamples, 1}));

        /* step 2. repartition in target space */
        final JavaPairRDD<LinearlySpacedIndexBlock, ImmutablePair<INDArray, INDArray>>
                blockifiedCopyRatioPosteriorResultsPairRDD = copyRatioPosteriorExpectationsPairRDD
                /* flat map to [target block, [sample index, [mean log copy ratio, var log copy ratio]]] */
                .flatMapToPair(dat -> targetBlocks.stream()
                        .map(tb -> new Tuple2<>(tb, new Tuple2<>(dat._1, ImmutablePair.of(
                                dat._2.getLogCopyRatioMeans(tb),
                                dat._2.getLogCopyRatioVariances(tb)))))
                        .iterator())
                /* combine the 1-to-many map from the previous step to a 1-to-1 map from target block to the list of
                 * sample-resolved copy ratio posteriors */
                .combineByKey(
                        /* recipe to create an singleton list */
                        Collections::singletonList,
                        /* recipe to add an element to the list */
                        (list, element) -> Stream.concat(list.stream(), Stream.of(element))
                                .collect(Collectors.toList()),
                        /* recipe to concatenate two lists */
                        (list1, list2) -> Stream.concat(list1.stream(), list2.stream()).collect(Collectors.toList()),
                        /* repartition with respect to target space blocks */
                        new HashPartitioner(numTargetBlocks))
                /* combine sample-resolved copy ratio posteriors in each target space block */
                .mapValues(list -> list.stream()
                        /* sort by sample index */
                        .sorted(Comparator.comparingInt(t -> t._1))
                        /* remove sample label */
                        .map(p -> p._2)
                        /* convert double[] to INDArray */
                        .map(t -> ImmutablePair.of(Nd4j.create(t.left), Nd4j.create(t.right)))
                        /* collect to a list */
                        .collect(Collectors.toList()))
                .mapValues(CoverageModelEMWorkspace::stackCopyRatioPosteriorDataForAllSamples);

        /* we do not need copy ratio expectations anymore */
        copyRatioPosteriorExpectationsPairRDD.unpersist();

        /* step 3. merge with computeRDD and update */
        computeRDD = computeRDD.join(blockifiedCopyRatioPosteriorResultsPairRDD)
                .mapValues(t -> t._1.cloneWithUpdatedCopyRatioPosteriors(t._2.left, t._2.right, admixingRatio));
        cacheWorkers("after E-step for copy ratio update");

        /* collect subroutine signals */
        final List<SubroutineSignal> sigs = mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal);

        final double errorNormInfinity = Collections.max(sigs.stream()
                .map(sig -> sig.<Double>get(StandardSubroutineSignals.RESIDUAL_ERROR_NORM))
                .collect(Collectors.toList()));

        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .build();
    }

    /**
     * Creates a {@link JavaPairRDD} of (sample index, emission data list)
     * @return a {@link JavaPairRDD}
     */
    private JavaPairRDD<Integer, List<CoverageModelCopyRatioEmissionData>> fetchCopyRatioEmissionDataSpark() {
        final int numSamples = this.numSamples;

        return computeRDD
                /* flat map workers a list of [sample index, [target block, emission data on target block]] */
                .flatMapToPair(tuple -> {
                    final LinearlySpacedIndexBlock tb = tuple._1;
                    final CoverageModelEMComputeBlock cb = tuple._2;
                    final List<List<CoverageModelCopyRatioEmissionData>> el = cb.getSampleCopyRatioLatentPosteriorData();
                    return IntStream.range(0, numSamples)
                            .mapToObj(si -> new Tuple2<>(si, new Tuple2<>(tb, el.get(si))))
                            .iterator();
                })
                /* combine elements with the same sample index */
                .combineByKey(
                        /* create a new list */
                        Collections::singletonList,
                        /* recipe to add an element to the list */
                        (list, element) -> Stream.concat(list.stream(), Stream.of(element))
                                .collect(Collectors.toList()),
                        /* recipe to concatenate two lists */
                        (list1, list2) -> Stream.concat(list1.stream(), list2.stream()).collect(Collectors.toList()),
                        /* repartition with respect to sample indices */
                        new HashPartitioner(numSamples))
                /* sort the [target block, emission data on target block] chunks for each sample into a single list */
                .mapValues(emissionBlocksList -> emissionBlocksList.stream() /* for each partition ... */
                        /* sort the data blocks */
                        .sorted(Comparator.comparingInt(t -> t._1.getBegIndex()))
                        /* remove the LinearlySpacedIndexBlock keys from the sorted emissionBlocksList */
                        .map(p -> p._2)
                        /* flatten */
                        .flatMap(List::stream)
                        /* collect as a single list */
                        .collect(Collectors.toList()));
    }

    /**
     * Stack a list of (1 x T INDArray, 1 x T INDArray) pairs along the 0th axis of each INDArray
     * and returns a (STATE x T INDArray, STATE x T INDArray) pair
     *
     * @param perSampleData a list of (1 x T INDArray, 1 x T INDArray) pairs
     * @return a (STATE x T INDArray, STATE x T INDArray) pair
     */
    private static ImmutablePair<INDArray, INDArray> stackCopyRatioPosteriorDataForAllSamples(
            final List<ImmutablePair<INDArray, INDArray>> perSampleData) {
        return ImmutablePair.of(Nd4j.vstack(perSampleData.stream().map(p -> p.left).collect(Collectors.toList())),
                Nd4j.vstack(perSampleData.stream().map(p -> p.right).collect(Collectors.toList())));
    }

    /**
     * Fetch copy ratio segments from compute blocks (Spark implementation)
     *
     * @return a list of {@link CopyRatioHMMResults}
     */
    private List<List<HiddenStateSegmentRecord<STATE, Target>>> getCopyRatioSegmentsSpark() {
        /* local final member variables for lambda capture */
        final List<Target> processedTargetList = new ArrayList<>();
        processedTargetList.addAll(this.processedTargetList);
        final List<SexGenotypeData> processedSampleSexGenotypeData = new ArrayList<>();
        processedSampleSexGenotypeData.addAll(this.processedSampleSexGenotypeData);
        final List<String> processedSampleNameList = new ArrayList<>();
        processedSampleNameList.addAll(this.processedSampleNameList);
        final INDArray sampleReadDepths = Transforms.exp(sampleMeanLogReadDepths, true);
        final CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, STATE> copyRatioExpectationsCalculator =
                this.copyRatioExpectationsCalculator;
        final BiFunction<SexGenotypeData, Target, STATE> referenceStateFactory = this.referenceStateFactory;

        return fetchCopyRatioEmissionDataSpark()
                /* let the workers run fb and Viterbi */
                .mapPartitionsToPair(it -> {
                    final List<Tuple2<Integer, CopyRatioHMMResults<
                            CoverageModelCopyRatioEmissionData, STATE>>> newPartitionData = new ArrayList<>();
                    while (it.hasNext()) {
                        final Tuple2<Integer, List<CoverageModelCopyRatioEmissionData>> prevDatum = it.next();
                        final int sampleIndex = prevDatum._1;
                        final CopyRatioCallingMetadata copyRatioCallingMetadata = CopyRatioCallingMetadata.builder()
                                .sampleName(processedSampleNameList.get(sampleIndex))
                                .sampleSexGenotypeData(processedSampleSexGenotypeData.get(sampleIndex))
                                .sampleCoverageDepth(sampleReadDepths.getDouble(sampleIndex))
                                .emissionCalculationStrategy(EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                                .build();
                        newPartitionData.add(new Tuple2<>(sampleIndex,
                                copyRatioExpectationsCalculator.getCopyRatioHMMResults(copyRatioCallingMetadata,
                                        processedTargetList, prevDatum._2)));
                    }
                    return newPartitionData.iterator();
                }, true)
                /* run segmentation on each sample */
                .mapPartitionsToPair(it -> {
                    final List<Tuple2<Integer, List<HiddenStateSegmentRecord<STATE, Target>>>> newPartitionData = new ArrayList<>();
                    while (it.hasNext()) {
                        final Tuple2<Integer, CopyRatioHMMResults<CoverageModelCopyRatioEmissionData, STATE>>
                                prevDatum = it.next();
                        final int sampleIndex = prevDatum._1;
                        final CopyRatioHMMResults<CoverageModelCopyRatioEmissionData, STATE> result = prevDatum._2;
                        final HMMSegmentProcessor<CoverageModelCopyRatioEmissionData, STATE, Target> processor =
                                new HMMSegmentProcessor<>(
                                        Collections.singletonList(result.getMetaData().getSampleName()),
                                        Collections.singletonList(result.getMetaData().getSampleSexGenotypeData()),
                                        referenceStateFactory,
                                        Collections.singletonList(new HashedListTargetCollection<>(processedTargetList)),
                                        Collections.singletonList(result.getForwardBackwardResult()),
                                        Collections.singletonList(result.getViterbiResult()));
                        newPartitionData.add(new Tuple2<>(sampleIndex, processor.getSegmentsAsList()));
                    }
                    return newPartitionData.iterator();
                })
                .collect()
                .stream()
                /* sort by sample index */
                .sorted(Comparator.comparingInt(t -> t._1))
                /* get rid of sample index */
                .map(t -> t._2)
                /* collect to a list */
                .collect(Collectors.toList());
    }


    /**
     * M-step update of mean log bias ($m_t$)
     *
     * @return a {@link SubroutineSignal} object containing the update size
     */
    @UpdatesRDD @CachesRDD
    public SubroutineSignal updateMeanLogBias(final boolean neglectBiasCovariates) {
        mapWorkers(cb -> cb
                .cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.M_STEP_M)
                .cloneWithUpdatedMeanLogBias(neglectBiasCovariates));
        cacheWorkers("after M-step for target mean bias");
        /* accumulate error from all nodes */
        final double errorNormInfinity = Collections.max(
                mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal).stream()
                        .map(sig -> sig.<Double>get(StandardSubroutineSignals.RESIDUAL_ERROR_NORM))
                        .collect(Collectors.toList()));

        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .build();
    }

    /**
     * M-step for target unexplained variance ($\Psi_t$)
     *
     * @return a {@link SubroutineSignal} object containing information about the solution
     */
    @UpdatesRDD @CachesRDD
    public SubroutineSignal updateTargetUnexplainedVariance(@Nonnull final CoverageModelArgumentCollection.TargetSpecificVarianceUpdateMode psiUpdateMode) {
        final int psiMaxIterations = config.getTargetSpecificVarianceMaxIterations();
        final double psiAbsoluteTolerance = config.getTargetSpecificVarianceAbsoluteTolerance();
        final double psiRelativeTolerance = config.getTargetSpecificVarianceRelativeTolerance();
        final double psiUpperLimit = config.getTargetSpecificVarianceUpperLimit();
        final int psiSolverNumBisections = config.getTargetSpecificVarianceSolverNumBisections();
        final int psiSolverRefinementDepth = config.getTargetSpecificVarianceSolverRefinementDepth();
        final int psiSolverNumThreads = config.getTargetSpecificVarianceSolverNumThreads();

        logger.debug("Target-specific variance update mode: " + config.getTargetSpecificVarianceUpdateMode().name());
        switch (psiUpdateMode) {
            case TARGET_RESOLVED: /* done on the compute blocks */
                mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.M_STEP_PSI)
                        .cloneWithUpdatedTargetUnexplainedVarianceTargetResolved(psiMaxIterations, psiUpperLimit, psiAbsoluteTolerance,
                                psiRelativeTolerance, psiSolverNumBisections, psiSolverRefinementDepth, psiSolverNumThreads));
                break;

            case ISOTROPIC: /* done on the driver node */
                return updateTargetUnexplainedVarianceIsotropic();

            default:
                throw new RuntimeException("Illegal Psi solver type.");
        }

        cacheWorkers("after M-step for target unexplained variance");

        /* accumulate error from all workers */
        final List<SubroutineSignal> signalList = mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal);
        final double errorNormInfinity = Collections.max(signalList.stream()
                .map(sig -> sig.<Double>get(StandardSubroutineSignals.RESIDUAL_ERROR_NORM))
                .collect(Collectors.toList()));
        final int maxIterations = Collections.max(signalList.stream()
                .map(sig -> sig.<Integer>get(StandardSubroutineSignals.ITERATIONS))
                .collect(Collectors.toList()));
        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .put(StandardSubroutineSignals.ITERATIONS, maxIterations) /* for uniformity */
                .build();
    }

    /**
     * M-step update of unexplained variance in the isotropic mode
     *
     * @return a {@link SubroutineSignal} object containing "error_norm" and "iterations" fields
     */
    @UpdatesRDD @CachesRDD
    private SubroutineSignal updateTargetUnexplainedVarianceIsotropic() {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.M_STEP_PSI));
        cacheWorkers("after M-step update of isotropic unexplained variance initialization");

        final double oldIsotropicTargetSpecificVariance = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t, 1)
                .meanNumber().doubleValue();

        final UnivariateFunction objFunc = psi -> mapWorkersAndReduce(cb ->
                cb.calculateSampleTargetSummedTargetSpecificVarianceObjectiveFunction(psi), (a, b) -> a + b);
        final UnivariateFunction meritFunc = psi -> mapWorkersAndReduce(cb ->
                cb.calculateSampleTargetSummedTargetSpecificVarianceMeritFunction(psi), (a, b) -> a + b);

        final RobustBrentSolver solver = new RobustBrentSolver(config.getTargetSpecificVarianceRelativeTolerance(),
                config.getTargetSpecificVarianceAbsoluteTolerance(), CoverageModelGlobalConstants.DEFAULT_FUNCTION_EVALUATION_ACCURACY,
                meritFunc, config.getTargetSpecificVarianceSolverNumBisections(), config.getTargetSpecificVarianceSolverRefinementDepth());
        double newIsotropicTargetSpecificVariance;
        try {
            newIsotropicTargetSpecificVariance = solver.solve(config.getTargetSpecificVarianceMaxIterations(), objFunc, 0, config.getTargetSpecificVarianceUpperLimit());
        } catch (NoBracketingException e) {
            logger.warn("Root of M-step optimality equation for isotropic unexplained variance could be bracketed");
            newIsotropicTargetSpecificVariance = oldIsotropicTargetSpecificVariance;
        } catch (TooManyEvaluationsException e) {
            logger.warn("Too many evaluations -- increase the number of root-finding iterations for the M-step update" +
                    " of unexplained variance");
            newIsotropicTargetSpecificVariance = oldIsotropicTargetSpecificVariance;
        }

        /* update the compute block(s) */
        final double errNormInfinity = FastMath.abs(newIsotropicTargetSpecificVariance - oldIsotropicTargetSpecificVariance);
        final int maxIterations = solver.getEvaluations();
        final double finalizedNewIsotropicTargetSpecificVariance = newIsotropicTargetSpecificVariance;
        mapWorkers(cb -> cb.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t,
                Nd4j.zeros(1, cb.getTargetSpaceBlock().getNumElements()).addi(finalizedNewIsotropicTargetSpecificVariance)));
        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errNormInfinity)
                .put(StandardSubroutineSignals.ITERATIONS, maxIterations).build();
    }

    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateBiasCovariatesARDCoefficients() {
        if (!ardEnabled) {
            throw new IllegalStateException("ARD update requested but it is disabled.");
        }
        final double ARD_ADMIXING_RATIO = 1.0;

        final INDArray alpha_l_inv_new = mapWorkersAndReduce(
                CoverageModelEMComputeBlock::getBiasCovariatesSecondMomentPosteriorsPartialTargetSum,
                INDArray::add);
        final double[] alpha_l_new_array = Nd4j.zeros(numLatents).addi(numTargets)
                .divi(alpha_l_inv_new).data().asDouble();
        final double[] alpha_l_old_array = biasCovariatesARDCoefficients.data().asDouble();
        final double[] alpha_l_new_admixed_array = IntStream.range(0, numLatents)
                .mapToDouble(li -> ARD_ADMIXING_RATIO * FastMath.min(FastMath.abs(alpha_l_new_array[li]),
                        config.getMaxARDPrecision()) + (1 - ARD_ADMIXING_RATIO) * alpha_l_old_array[li])
                .toArray();

        /* update */
        final INDArray alpha_l_new = Nd4j.create(alpha_l_new_admixed_array, new int[] {1, numLatents});

        /* info log to stdout */
        logger.info("Log ARD coefficients: " + Transforms.log(alpha_l_new, true).toString());

        final double errorNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                Transforms.log(alpha_l_new, true).subi(Transforms.log(biasCovariatesARDCoefficients, true)));

        /* update the worker copy */
        mapWorkers(cb -> cb.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.alpha_l,
                alpha_l_new));

        /* update the driver copy */
        biasCovariatesARDCoefficients.assign(alpha_l_new);

        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .build();
    }

    /**
     * E-step update of bias covariates
     *
     * @return a {@link SubroutineSignal} object containing "error_norm"
     */
    @UpdatesRDD @CachesRDD
    public SubroutineSignal updateBiasCovariates(final double admixingRatio) {
        /* perform the M-step update */
        final SubroutineSignal sig;
        if (!config.fourierRegularizationEnabled()) {
            sig = updateBiasCovariatesUnregularized(admixingRatio);
        } else {
            sig = updateBiasCovariatesRegularized(admixingRatio);
        }
        return sig;
    }

    @UpdatesRDD @CachesRDD
    public SubroutineSignal updateBiasCovariates() {
        return updateBiasCovariates(config.getMeanFieldAdmixingRatio());
    }

    /**
     * E-step update of bias covariates w/o regularization
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    @UpdatesRDD @CachesRDD
    private SubroutineSignal updateBiasCovariatesUnregularized(final double admixingRatio) {
        mapWorkers(cb -> cb
                .cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_W_UNREG)
                .cloneWithUpdatedBiasCovariatesUnregularized(admixingRatio));
        cacheWorkers("after E-step update of bias covariates w/o regularization");

        /* accumulate error from all nodes */
        final double errorNormInfinity = Collections.max(
                mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal).stream()
                        .map(sig -> sig.<Double>get(StandardSubroutineSignals.RESIDUAL_ERROR_NORM))
                        .collect(Collectors.toList()));

        return SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .build();
    }

    /**
     * E-step update of bias covariates w/ regularization (local implementation)
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    @UpdatesRDD @EvaluatesRDD @CachesRDD
    private SubroutineSignal updateBiasCovariatesRegularized(final double admixingRatio) {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_W_REG));
        cacheWorkers("after E-step update of bias covariates w/ regularization");

        final INDArray W_tl_old = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl, 0);
        final INDArray v_tl = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.v_tl, 0);

        /* initialize the linear operators */
        final GeneralLinearOperator<INDArray> linop, precond;
        final ImmutablePair<GeneralLinearOperator<INDArray>, GeneralLinearOperator<INDArray>> ops =
                getBiasCovariatesRegularizedLinearOperators();
        linop = ops.left;
        precond = ops.right;

        /* initialize the iterative solver */
        final IterativeLinearSolverNDArray iterSolver = new IterativeLinearSolverNDArray(linop, v_tl, precond,
                config.getWAbsoluteTolerance(), config.getWRelativeTolerance(), config.getWMaxIterations(),
                x -> x.normmaxNumber().doubleValue(), /* norm */
                (x, y) -> x.mul(y).sumNumber().doubleValue(), /* inner product */
                true);

        /* solve */
        long startTime = System.nanoTime();
        final SubroutineSignal sig = iterSolver.solveUsingPreconditionedConjugateGradient(W_tl_old);
        linop.cleanupAfter();
        precond.cleanupAfter();
        long endTime = System.nanoTime();
        logger.debug("CG execution time for solving the regularized M-step update equation for bias covariates" +
                (double)(endTime - startTime)/1000000 + " ms");

        /* check the exit status of the solver and push the new W to workers */
        final ExitStatus exitStatus = sig.get(StandardSubroutineSignals.EXIT_STATUS);
        if (exitStatus == ExitStatus.FAIL_MAX_ITERS) {
            logger.warn("CG iterations for M-step update of bias covariates did not converge. Increase maximum iterations" +
                    " and/or decrease absolute/relative error tolerances");
        }
        final int iters = sig.<Integer>get(StandardSubroutineSignals.ITERATIONS);
        final INDArray W_tl_new = sig.get(StandardSubroutineSignals.SOLUTION);

        switch (config.getBiasCovariatesComputeNodeCommunicationPolicy()) {
            case BROADCAST_HASH_JOIN:
                pushToWorkers(mapINDArrayToBlocks(W_tl_new), (W, cb) ->
                        cb.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl,
                                W.get(cb.getTargetSpaceBlock())));
                break;

            case RDD_JOIN:
                joinWithWorkersAndMap(chopINDArrayToBlocks(W_tl_new),
                        p -> p._1.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl,
                                p._2));
                break;

            default:
                throw new GATKException.ShouldNeverReachHereException("Unknown communication policy for M-step update" +
                        " of bias covariates");
        }

        /* update F[W] */
        updateFilteredBiasCovariates(W_tl_new);

        final double errorNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(W_tl_new.sub(W_tl_old));

        /* send the signal to workers for consistency */
        final SubroutineSignal newSig = SubroutineSignal.builder()
                .put(StandardSubroutineSignals.RESIDUAL_ERROR_NORM, errorNormInfinity)
                .put(StandardSubroutineSignals.ITERATIONS, iters)
                .build();
        mapWorkers(cb -> cb.cloneWithUpdatedSignal(newSig));
        return newSig;
    }

    /**
     * Returns the pair of (linear operator, preconditioner) invoked in the M-step update of bias
     * covariates in w/ regularization
     *
     * @return a pair of linear operators
     */
    @EvaluatesRDD
    private ImmutablePair<GeneralLinearOperator<INDArray>, GeneralLinearOperator<INDArray>> getBiasCovariatesRegularizedLinearOperators() {
        final INDArray Q_ll, Q_tll, Z_ll;
        final GeneralLinearOperator<INDArray> linop, precond;
        final FourierLinearOperatorNDArray regularizerFourierLinearOperator = createRegularizerFourierLinearOperator();

        switch (config.getWSolverType()) {
            case LOCAL:
                /* fetch the required INDArrays */
                Q_ll = mapWorkersAndReduce(cb -> cb.getINDArrayFromCache(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.sum_Q_ll),
                        INDArray::add).div(numTargets);
                Q_tll = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Q_tll, 0);
                Z_ll = sampleBiasLatentPosteriorSecondMoments.sum(0);

                /* instantiate the local implementation of linear operators */
                linop = new CoverageModelWLinearOperatorLocal(Q_tll, Z_ll, regularizerFourierLinearOperator);
                precond = new CoverageModelWPreconditionerLocal(Q_ll, Z_ll, regularizerFourierLinearOperator, numTargets);

                return ImmutablePair.of(linop, precond);

            case SPARK:
                if (!sparkContextIsAvailable) {
                    throw new UserException("The Spark W solver is only available in the Spark mode");
                }
                /* fetch the required INDArrays */
                Q_ll = mapWorkersAndReduce(cb -> cb.getINDArrayFromCache(CoverageModelEMComputeBlock
                        .CoverageModelICGCacheNode.sum_Q_ll), INDArray::add).div(numTargets);
                Z_ll = sampleBiasLatentPosteriorSecondMoments.sum(0);

                /* instantiate the spark implementation of linear operators */
                linop = new CoverageModelWLinearOperatorSpark(Z_ll, regularizerFourierLinearOperator,
                        numTargets, ctx, computeRDD, targetBlocks);
                precond = new CoverageModelWPreconditionerSpark(Q_ll, Z_ll, regularizerFourierLinearOperator,
                        numTargets, ctx, numTargetBlocks);

                return ImmutablePair.of(linop, precond);

            default:
                throw new IllegalArgumentException("The solver type for M-step update of bias covariates" +
                        " is not properly set");
        }
    }

    /**
     * Fetch the log likelihood from compute block(s)
     *
     * Include the ARD coefficients
     *
     * @return log likelihood normalized per sample per target
     */
    @EvaluatesRDD
    public double getLogLikelihood() {
        /* this is normalized by active targets of each sample */
        final double logLikelihoodEmissionAndChainPerTarget = Arrays.stream(getLogLikelihoodPerSample()).reduce((a, b) -> a + b)
                .orElse(Double.NaN);

        /* contribution from ARD */
        final double logLikelihoodARDPerTarget;
        if (ardEnabled && config.includeARDInLogLikelihood()) {
            final INDArray biasCovariatesSecondMoment = mapWorkersAndReduce(
                    CoverageModelEMComputeBlock::getBiasCovariatesSecondMomentPosteriorsPartialTargetSum, INDArray::add);
            logLikelihoodARDPerTarget = 0.5 * Transforms.log(biasCovariatesARDCoefficients, true).sumNumber().doubleValue()
                    - 0.5 * FastMath.log(2 * FastMath.PI)
                    - 0.5 * biasCovariatesARDCoefficients.mul(biasCovariatesSecondMoment).sumNumber().doubleValue() / numTargets;
        } else {
            logLikelihoodARDPerTarget = 0;
        }

        /* the final result is normalized for number of samples */
        double logLikelihoodTotal = (logLikelihoodEmissionAndChainPerTarget + logLikelihoodARDPerTarget) / numSamples;

        /* log history */
        logLikelihoodHistory.add(logLikelihoodTotal);
        if (ardEnabled) {
            biasCovariatesARDCoefficientsHistory.add(biasCovariatesARDCoefficients.dup());
        }

        return logLikelihoodTotal;
    }

    /**
     * Fetch the log likelihood from compute block(s)
     *
     * @return log likelihood normalized per target
     */
    @EvaluatesRDD @CachesRDD
    public double[] getLogLikelihoodPerSample() {
        updateLogLikelihoodCaches();

        /* contribution from latent bias prior */
        final INDArray biasPriorContrib = Nd4j.zeros(new int[] {numSamples, 1});
        if (biasCovariatesEnabled) {
            for (int li = 0; li < numLatents; li++) {
                biasPriorContrib.addi(sampleBiasLatentPosteriorSecondMoments.get(NDArrayIndex.all(),
                        NDArrayIndex.point(li), NDArrayIndex.point(li)));
            }
            biasPriorContrib.addi(0.5 * numLatents * FastMath.log(2 * FastMath.PI)).muli(-0.5);
        }

        /* contribution from emission edges */
        final CoverageModelEMComputeBlock.CoverageModelICGCacheNode logLikelihoodKey;
        if (config.fourierRegularizationEnabled()) {
            logLikelihoodKey = CoverageModelEMComputeBlock.CoverageModelICGCacheNode.loglike_reg;
        } else {
            logLikelihoodKey = CoverageModelEMComputeBlock.CoverageModelICGCacheNode.loglike_unreg;
        }
        final INDArray emissionContrib = mapWorkersAndReduce(cb -> cb.getINDArrayFromCache(logLikelihoodKey),
                INDArray::add);

        /* number of unmasked (active for learning) targets per sample */
        final INDArray sum_M = mapWorkersAndReduce(cb -> cb.getINDArrayFromCache(
                CoverageModelEMComputeBlock.CoverageModelICGCacheNode.sum_M_s), INDArray::add);

        final INDArray logLikelihoodEmissionPerTarget = emissionContrib.add(biasPriorContrib).divi(sum_M);
        final INDArray logLikelihoodChainPerTarget = sampleLogChainPosteriors.div(numTargets);

        return logLikelihoodEmissionPerTarget.addi(logLikelihoodChainPerTarget).data().asDouble();
    }

    /**
     * Updates log likelihood caches on compute blocks
     */
    @UpdatesRDD @CachesRDD
    public void updateLogLikelihoodCaches() {
        if (config.fourierRegularizationEnabled()) {
            mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock
                    .CoverageModelICGCacheTag.LOGLIKE_REG));
        } else {
            mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock
                    .CoverageModelICGCacheTag.LOGLIKE_UNREG));
        }
        cacheWorkers("after updating log likelihood caches");
    }

    /**
     * Instantiate compute block(s). If Spark is disabled, a single {@link CoverageModelEMComputeBlock} is
     * instantiated. Otherwise, a {@link JavaPairRDD} of compute nodes will be created.
     */
    private void instantiateWorkers() {
        if (sparkContextIsAvailable) {
            /* initialize the RDD */
            logger.info("Initializing an RDD of compute blocks");
            computeRDD = ctx.parallelizePairs(targetBlockStream()
                    .map(tb -> new Tuple2<>(tb, new CoverageModelEMComputeBlock(tb, numSamples, numLatents, ardEnabled)))
                    .collect(Collectors.toList()), numTargetBlocks)
                    .partitionBy(new HashPartitioner(numTargetBlocks))
                    .cache();
        } else {
            logger.info("Initializing a local compute block");
            localComputeBlock = new CoverageModelEMComputeBlock(targetBlocks.get(0), numSamples, numLatents, ardEnabled);
        }
        prevCheckpointedComputeRDD = null;
        cacheCallCounter = 0;
    }

    /**
     * A generic function for handling a blockified list of objects to their corresponding compute nodes
     *
     * If Spark is enabled:
     *
     *      Joins an instance of {@code List<Tuple2<LinearlySpacedIndexBlock, V>>} with {@link #computeRDD}, calls the provided
     *      map {@code mapper} on the RDD, and the reference to the old RDD will be replaced with the new RDD.
     *
     * If Spark is disabled:
     *
     *      Only a single target-space block is assumed, such that {@code data} is a singleton. The map function
     *      {@code mapper} will be called on the value contained in {@code data} and {@link #localComputeBlock}, and
     *      the old instance of {@link CoverageModelEMComputeBlock} is replaced with the new instance returned
     *      by {@code mapper.}
     *
     * @param data the list to joined and mapped together with the compute block(s)
     * @param mapper a mapper binary function that takes a compute block together with an object of type {@code V} and
     *               returns a new compute block
     * @param <V> the type of the object to the broadcasted
     */
    @UpdatesRDD
    private <V> void joinWithWorkersAndMap(@Nonnull final List<Tuple2<LinearlySpacedIndexBlock, V>> data,
                                           @Nonnull final Function<Tuple2<CoverageModelEMComputeBlock, V>, CoverageModelEMComputeBlock> mapper) {
        if (sparkContextIsAvailable) {
            final JavaPairRDD<LinearlySpacedIndexBlock, V> newRDD =
                    ctx.parallelizePairs(data, numTargetBlocks).partitionBy(new HashPartitioner(numTargetBlocks));
            computeRDD = computeRDD.join(newRDD).mapValues(mapper);
        } else {
            try {
                Utils.validateArg(data.size() == 1, "Only a single data block is expected in the local mode");
                localComputeBlock = mapper.call(new Tuple2<>(localComputeBlock, data.get(0)._2));
            } catch (Exception e) {
                throw new RuntimeException("Can not apply the map function to the local compute block: " + e.getMessage());
            }
        }
    }

    /**
     * Calls a map function on the compute block(s) and returns the new compute block(s)
     *
     * If Spark is enabled:
     *
     *      The map is applied on the values of {@link #computeRDD}, the reference to the old RDD will be replaced
     *      by the new RDD, and original partitioning is retained
     *
     * If Spark is disabled:
     *
     *      Only a single target-space block is assumed; the map is applied to {@link #localComputeBlock} and the
     *      reference is updated accordingly
     *
     * @param mapper a map from {@link CoverageModelEMComputeBlock} onto itself
     */
    @UpdatesRDD
    private void mapWorkers(@Nonnull final Function<CoverageModelEMComputeBlock, CoverageModelEMComputeBlock> mapper) {
        if (sparkContextIsAvailable) {
            computeRDD = computeRDD.mapValues(mapper);
        } else {
            try {
                localComputeBlock = mapper.call(localComputeBlock);
            } catch (final Exception ex) {
                throw new RuntimeException("Can not apply the map function to the local compute block", ex);
            }
        }
    }

    /**
     * Calls a generic map function on compute block(s) and collects the values to a {@link List}
     *
     * If Spark is enabled:
     *
     *      The size of the list is the same as the number of elements in the RDD
     *
     * If Spark is disabled:
     *
     *      The list will be singleton
     *
     * @param mapper a map function from {@link CoverageModelEMComputeBlock} to a generic type
     * @param <V> the return type of the map function
     * @return a list of collected mapped values
     */
    @EvaluatesRDD
    private <V> List<V> mapWorkersAndCollect(@Nonnull final Function<CoverageModelEMComputeBlock, V> mapper) {
        if (sparkContextIsAvailable) {
            return computeRDD.values().map(mapper).collect();
        } else {
            try {
                return Collections.singletonList(mapper.call(localComputeBlock));
            } catch (final Exception ex) {
                throw new RuntimeException("Can not apply the map function to the local compute block", ex);
            }
        }
    }

    /**
     * A generic map-reduce step on the compute block(s)
     *
     * If Spark is enabled:
     *
     *      Map the RDD by {@code mapper} and reduce by {@code reducer}
     *
     * If Spark is disabled:
     *
     *      Call on the map and reduce on {@link #localComputeBlock}
     *
     * @param mapper a map from {@link CoverageModelEMComputeBlock} to a generic type
     * @param reducer a generic symmetric reducer binary function from (V, V) -> V
     * @param <V> the type of the reduction
     * @return the result of map-reduce
     */
    @EvaluatesRDD
    private <V> V mapWorkersAndReduce(@Nonnull final Function<CoverageModelEMComputeBlock, V> mapper,
                                      @Nonnull final Function2<V, V, V> reducer) {
        if (sparkContextIsAvailable) {
            return computeRDD.values().map(mapper).reduce(reducer);
        } else {
            try {
                return mapper.call(localComputeBlock);
            } catch (final Exception ex) {
                throw new RuntimeException("Can not apply the map function to the local compute block", ex);
            }

        }
    }

    /**
     * A generic function for broadcasting an object to all compute blocks
     *
     * If Spark is enabled:
     *
     *      A {@link Broadcast} will be created from {@param obj} and will be "received" by the compute nodes by calling
     *      {@param pusher}. A reference to the updated RDD will replace the old RDD.
     *
     * If Spark is disabled:
     *
     *      The {@param pusher} function will be called together with {@param obj} and {@link #localComputeBlock}
     *
     * @param obj te object to broadcast
     * @param pusher a map from (V, {@link CoverageModelEMComputeBlock}) -> {@link CoverageModelEMComputeBlock} that
     *               updates the compute block with the broadcasted value
     * @param <V> the type of the broadcasted object
     */
    @UpdatesRDD
    private <V> void pushToWorkers(@Nonnull final V obj,
                                   @Nonnull final Function2<V, CoverageModelEMComputeBlock, CoverageModelEMComputeBlock> pusher) {
        if (sparkContextIsAvailable) {
            final Broadcast<V> broadcastedObj = ctx.broadcast(obj);
            final Function<CoverageModelEMComputeBlock, CoverageModelEMComputeBlock> mapper =
                    cb -> pusher.call(broadcastedObj.value(), cb);
            mapWorkers(mapper);
        } else {
            try {
                localComputeBlock = pusher.call(obj, localComputeBlock);
            } catch (final Exception ex) {
                throw new RuntimeException("Can not apply the map function to the local compute block", ex);
            }
        }
    }

    /**
     * If Spark is enabled, caches the RDD of compute block(s). Otherwise, it does nothing.
     *
     * @param where a message provided by the method that calls this function
     */
    @CachesRDD
    public void cacheWorkers(final String where) {
        if (sparkContextIsAvailable) {
            logger.debug("RDD caching requested (" + where + ")");
            computeRDD.persist(StorageLevel.MEMORY_AND_DISK_SER());
            cacheCallCounter++;
            if (!prevCachedComputeRDDDeque.isEmpty()) {
                prevCachedComputeRDDDeque.removeFirst().unpersist(true);
                prevCachedComputeRDDDeque.addLast(computeRDD);
            }
            if (config.isRDDCheckpointingEnabled()) {
                if (cacheCallCounter == config.getRDDCheckpointingInterval()) {
                    logger.debug("Checkpointing compute RDD...");
                    computeRDD.checkpoint();
                    if (prevCheckpointedComputeRDD != null) {
                        prevCheckpointedComputeRDD.unpersist(true);
                        prevCheckpointedComputeRDD = computeRDD;
                    }
                    cacheCallCounter = 0;
                }
            }
        }
    }

    /**
     * Fetches the blocks of a target-distributed {@link INDArray} of shape ({@link #numTargets}, ...)
     * and assembles them together by concatenating along {@param axis} (if spark is enabled)
     *
     * If Spark is disabled, it just fetches the {@link INDArray} from {@link #localComputeBlock}.
     *
     * @param key key of the array
     * @param axis axis to stack along
     * @return assembled array
     */
    @EvaluatesRDD @VisibleForTesting
    private INDArray fetchFromWorkers(final CoverageModelEMComputeBlock.CoverageModelICGCacheNode key, final int axis) {
        if (sparkContextIsAvailable) {
            return CoverageModelSparkUtils.assembleINDArrayBlocksFromRDD(computeRDD.mapValues(cb -> cb.getINDArrayFromCache(key)), axis);
        } else {
            return localComputeBlock.getINDArrayFromCache(key);
        }
    }

    /**
     * Partitions an {@link INDArray} along its first dimension and makes a key-value {@link List} of the blocks
     *
     * @param arr the input array
     * @return list of key-value blocks
     */
    private List<Tuple2<LinearlySpacedIndexBlock, INDArray>> chopINDArrayToBlocks(final INDArray arr) {
        if (sparkContextIsAvailable) {
            return CoverageModelSparkUtils.partitionINDArrayToList(targetBlocks, arr);
        } else {
            return Collections.singletonList(new Tuple2<>(targetBlocks.get(0), arr));
        }
    }

    /**
     * Partitions an {@link INDArray} along its first dimension and makes a map
     *
     * @param arr the input array
     * @return list of key-value blocks
     */
    private Map<LinearlySpacedIndexBlock, INDArray> mapINDArrayToBlocks(final INDArray arr) {
        if (sparkContextIsAvailable) {
            return CoverageModelSparkUtils.partitionINDArrayToMap(targetBlocks, arr);
        } else {
            return Collections.singletonMap(targetBlocks.get(0), arr);
        }
    }

    /**
     * Partitions two INDArrays simultaneously into target-space blocks along the 0th axis
     *
     * @param arr1 first INDArray
     * @param arr2 second INDArray
     * @return a map from target-space blocks to partitioned pairs
     */
    private Map<LinearlySpacedIndexBlock, ImmutablePair<INDArray, INDArray>> mapINDArrayPairToBlocks(final INDArray arr1,
                                                                                                     final INDArray arr2) {
        if (sparkContextIsAvailable) {
            final Map<LinearlySpacedIndexBlock, INDArray> map1 =
                    CoverageModelSparkUtils.partitionINDArrayToMap(targetBlocks, arr1);
            final Map<LinearlySpacedIndexBlock, INDArray> map2 =
                    CoverageModelSparkUtils.partitionINDArrayToMap(targetBlocks, arr2);
            final Map<LinearlySpacedIndexBlock, ImmutablePair<INDArray, INDArray>> res = new HashMap<>();
            targetBlockStream().forEach(tb -> res.put(tb, ImmutablePair.of(map1.get(tb), map2.get(tb))));
            return res;
        } else {
            return Collections.singletonMap(targetBlocks.get(0), ImmutablePair.of(arr1, arr2));
        }
    }

    /**
     * Returns an {@link IntStream} of sample indices
     *
     * @return {@link IntStream}
     */
    private IntStream sampleIndexStream() { return IntStream.range(0, numSamples); }

    /**
     * Returns a {@link Stream< LinearlySpacedIndexBlock >} of target-space blocks
     *
     * @return {@link Stream< LinearlySpacedIndexBlock >}
     */
    private Stream<LinearlySpacedIndexBlock> targetBlockStream() { return targetBlocks.stream(); }

    /**
     * Calls gc() on all compute blocks
     *
     * Refer to te implementation docs of {@link CoverageModelEMComputeBlock#performGarbageCollection()}
     * for details.
     *
     */
    @EvaluatesRDD
    public void performGarbageCollection() {
        System.gc();
        if (sparkContextIsAvailable) {
            computeRDD.values().foreach(CoverageModelEMComputeBlock::performGarbageCollection);
        }
    }

    /**
     * Fetches mean log bias from compute blocks
     *
     * @return an {@link INDArray}
     */
    private INDArray fetchMeanLogBias() {
        return fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t, 1);
    }

    /**
     * Fetches total unexplained bias variance as a sample-target matrix from compute blocks
     *
     * @return an {@link INDArray}
     */
    private INDArray fetchTotalUnexplainedVariance() {
        return fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.tot_Psi_st, 1).transpose();
    }

    /**
     * Fetches total covariate bias [W.z]_{st} as a sample-target matrix from compute blocks
     *
     * @return an {@link INDArray}
     */
    private INDArray fetchTotalCovariateBiasPerSample() {
        return fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Wz_st, 1).transpose();
    }

    /**
     * Fetches target-specific unexplained bias variance from compute blocks
     *
     * @return an {@link INDArray}
     */
    private INDArray fetchTargetUnexplainedVariance() {
        return fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t, 1);
    }

    /**
     * Fetches mean bias covariates from compute blocks
     *
     * @return an {@link INDArray}
     */
    private INDArray fetchMeanBiasCovariates() {
        if (biasCovariatesEnabled) {
            return fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl, 0);
        } else {
            throw new IllegalStateException("Bias covariates are disabled but their means were fetched");
        }
    }

    /**
     * Fetches the maximum likelihood estimate of copy ratios and their precisions from compute blocks.
     * The result is output as a pair of target-by-sample matrices.
     *
     * @param logScale if true, the max likelihood estimate is reported in natural log scale
     * @return a pair of {@link INDArray}
     */
    private ImmutablePair<INDArray, INDArray> fetchCopyRatioMaxLikelihoodEstimateData(final boolean logScale) {
        final INDArray M_Psi_inv_st = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.M_Psi_inv_st, 1);
        final INDArray log_n_st = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.log_n_st, 1);
        final INDArray m_t = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t, 1);

        /* calculate the required quantities */
        final INDArray copyRatioMaxLikelihoodEstimate;
        if (biasCovariatesEnabled) {
            final INDArray Wz_st = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Wz_st, 1);
            copyRatioMaxLikelihoodEstimate = log_n_st.sub(Wz_st).subiRowVector(m_t).subiColumnVector(sampleMeanLogReadDepths);
        } else {
            copyRatioMaxLikelihoodEstimate = log_n_st.subRowVector(m_t).subiColumnVector(sampleMeanLogReadDepths);
        }

        if (!logScale) {
            Transforms.exp(copyRatioMaxLikelihoodEstimate, false);
        }

        return ImmutablePair.of(copyRatioMaxLikelihoodEstimate.transpose(), M_Psi_inv_st.transpose());
    }

    /**
     * Fetches the Viterbi copy ratio (or copy number) states as a target-sample matrix
     *
     * @return an {@link INDArray}
     */
    protected INDArray getViterbiAsNDArray(final List<List<HiddenStateSegmentRecord<STATE, Target>>> segments) {
        final INDArray res = Nd4j.create(numSamples, numTargets);
        final TargetCollection<Target> targetCollection = new HashedListTargetCollection<>(processedTargetList);
        for (int si = 0; si < numSamples; si++) {
            final SexGenotypeData sampleSexGenotype = processedSampleSexGenotypeData.get(si);
            /* start with all ref */
            final double[] calls = IntStream.range(0, numTargets).mapToDouble(
                    ti -> referenceStateFactory.apply(sampleSexGenotype, processedTargetList.get(ti)).getScalar())
                    .toArray();
            /* go through segments and mutate ref calls as necessary */
            segments.get(si).forEach(seg -> {
                final IndexRange range = targetCollection.indexRange(seg.getSegment());
                final double copyRatio = seg.getSegment().getCall().getScalar();
                for (int ti = range.from; ti < range.to; ti++) {
                    calls[ti] = copyRatio;
                }
            });
            res.getRow(si).assign(Nd4j.create(calls, new int[] {1, numTargets}));
        }
        return res.transpose();
    }

    /**
     * Saves the model to disk
     *
     * @param outputPath path to write to the model
     */
    public void writeModel(@Nonnull final String outputPath) {
        logger.info("Saving the model to disk...");
        CoverageModelParameters.write(new CoverageModelParameters(
                processedTargetList, fetchMeanLogBias(), fetchTargetUnexplainedVariance(),
                biasCovariatesEnabled ? fetchMeanBiasCovariates() : null,
                ardEnabled ? biasCovariatesARDCoefficients : null), outputPath);
    }

    /**
     * Saves posteriors to disk
     *
     * @param outputPath path to write posteriors
     * @param verbosityLevel verbosity level
     */
    public void writePosteriors(final String outputPath, final PosteriorVerbosityLevel verbosityLevel) {
        /* create output directory if it doesn't exist */
        createOutputPath(outputPath);

        writeTargets(outputPath);
        writeReadDepthPosteriors(outputPath);
        writeLogLikelihoodPosteriors(outputPath);
        writeSampleSpecificUnexplainedVariancePosteriors(outputPath);
        if (biasCovariatesEnabled) {
            writeBiasLatentPosteriors(outputPath);
        }
        if (ardEnabled) {
            writeBiasCovariatesARDHistory(outputPath);
        }
        writeCopyRatioPosteriors(outputPath);
        writeCopyRatioMaxLikelihoodEstimates(outputPath);

        if (verbosityLevel.equals(PosteriorVerbosityLevel.EXTENDED)) {
            writeExtendedPosteriors(outputPath);
        }
    }

    /**
     * Create output path if non-existent
     *
     * @param outputPath the output path
     */
    private void createOutputPath(final String outputPath) {
        final File outputPathFile = new File(outputPath);
        if (!outputPathFile.exists()) {
            if (!outputPathFile.mkdirs()) {
                throw new UserException.CouldNotCreateOutputFile(outputPathFile, "Could not create the output directory");
            }
        }
    }

    /**
     * Save targets to disk
     *
     * @param outputPath the output path
     */
    protected void writeTargets(final String outputPath) {
        logger.info("Saving targets...");
        final File targetListFile = new File(outputPath, CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE);
        TargetWriter.writeTargetsToFile(targetListFile, processedTargetList);
    }

    /**
     * Saves copy ratio (or copy number) max likelihood estimates to disk
     *
     * @param outputPath the output path
     */
    protected void writeCopyRatioMaxLikelihoodEstimates(final String outputPath) {
        final ImmutablePair<INDArray, INDArray> copyRatioMLEData = fetchCopyRatioMaxLikelihoodEstimateData(false);
        final List<String> sampleNames = processedReadCounts.columnNames();
        final List<String> targetNames = processedReadCounts.targets().stream()
                .map(Target::getName).collect(Collectors.toList());

        final File copyRatioMLEFile = new File(outputPath, CoverageModelGlobalConstants.COPY_RATIO_MAX_LIKELIHOOD_ESTIMATES_FILENAME);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(copyRatioMLEData.left, copyRatioMLEFile, "COPY_RATIO_MAX_LIKELIHOOD_ESTIMATES",
                targetNames, sampleNames);

        final File copyRatioPrecisionFile = new File(outputPath, CoverageModelGlobalConstants.COPY_RATIO_PRECISION_FILENAME);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(copyRatioMLEData.right, copyRatioPrecisionFile, "LOG_COPY_RATIO_PRECISIONS",
                targetNames, sampleNames);
    }

    /**
     * Saves read depth posteriors to disk
     *
     * @param outputPath the output path
     */
    protected void writeReadDepthPosteriors(final String outputPath) {
        logger.info("Saving read depth posteriors...");
        final List<String> sampleNames = processedReadCounts.columnNames();
        final INDArray combinedReadDepthPosteriors = Nd4j.hstack(sampleMeanLogReadDepths, sampleVarLogReadDepths);
        final File sampleReadDepthPosteriorsFile = new File(outputPath, CoverageModelGlobalConstants.SAMPLE_READ_DEPTH_POSTERIORS_FILENAME);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(combinedReadDepthPosteriors, sampleReadDepthPosteriorsFile,
                "SAMPLE_NAME", sampleNames, Arrays.asList("READ_DEPTH_MEAN", "READ_DEPTH_VAR"));
    }

    /**
     * Saves model mog likelihood posteriors to disk
     *
     * @param outputPath the output path
     */
    protected void writeLogLikelihoodPosteriors(final String outputPath) {
        logger.info("Saving sample log likelihoods...");
        final List<String> sampleNames = processedReadCounts.columnNames();
        final File sampleLogLikelihoodsFile = new File(outputPath, CoverageModelGlobalConstants.SAMPLE_LOG_LIKELIHOODS_FILENAME);
        final INDArray sampleLogLikelihoods = Nd4j.create(getLogLikelihoodPerSample(), new int[] {numSamples, 1});
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(sampleLogLikelihoods, sampleLogLikelihoodsFile,
                "SAMPLE_NAME", sampleNames, Collections.singletonList("LOG_LIKELIHOOD"));

        /* write log likelihood history */
        if (logLikelihoodHistory.size() > 0) {
            final File logLikelihoodHistoryFile = new File(outputPath, CoverageModelGlobalConstants.LOG_LIKELIHOODS_HISTORY_FILENAME);
            final INDArray logLikelihoodHistoryNDArray = Nd4j.create(logLikelihoodHistory.stream()
                    .mapToDouble(Double::valueOf).toArray(), new int[] {logLikelihoodHistory.size(), 1});
            Nd4jIOUtils.writeNDArrayMatrixToTextFile(logLikelihoodHistoryNDArray, logLikelihoodHistoryFile,
                    "LOG_LIKELIHOOD", null, null);
        }
    }

    /**
     * Saves sample-specific unexplained variance to disk
     *
     * @param outputPath the output path
     */
    protected void writeSampleSpecificUnexplainedVariancePosteriors(final String outputPath) {
        logger.info("Saving sample-specific unexplained variance posteriors...");
        final List<String> sampleNames = processedReadCounts.columnNames();
        final File sampleUnexplainedVarianceFile = new File(outputPath, CoverageModelGlobalConstants.SAMPLE_UNEXPLAINED_VARIANCE_FILENAME);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(sampleUnexplainedVariance, sampleUnexplainedVarianceFile,
                "SAMPLE_NAME", sampleNames, Collections.singletonList("SAMPLE_UNEXPLAINED_VARIANCE"));
    }

    /**
     * Saves bias latent posteriors E[z_{s\mu}] to disk
     *
     * @param outputPath the output path
     */
    protected void writeBiasLatentPosteriors(final String outputPath) {
        logger.info("Saving bias latent posteriors...");
        final List<String> sampleNames = processedReadCounts.columnNames();
        final File sampleBiasLatentPosteriorsFile = new File(outputPath,
                CoverageModelGlobalConstants.SAMPLE_BIAS_LATENT_POSTERIORS_FILENAME);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(sampleBiasLatentPosteriorFirstMoments, sampleBiasLatentPosteriorsFile,
                "SAMPLE_NAME", sampleNames, IntStream.range(0, numLatents)
                        .mapToObj(li -> String.format("PC_%d", li)).collect(Collectors.toList()));
    }

    protected void writeBiasCovariatesARDHistory(final String outputPath) {
        logger.info("Saving ARD history...");
        final File biasCovariatesARDHistoryFile = new File(outputPath,
                CoverageModelGlobalConstants.BIAS_COVARIATES_ARD_COEFFICIENTS_HISTORY_OUTPUT_FILE);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(Nd4j.vstack(biasCovariatesARDCoefficientsHistory),
                biasCovariatesARDHistoryFile, "ARD_HISTORY",
                IntStream.range(0, biasCovariatesARDCoefficientsHistory.size())
                        .mapToObj(i -> String.format("ITER_%d", i+1)).collect(Collectors.toList()),
                IntStream.range(0, numLatents)
                        .mapToObj(i -> String.format("PC_%d", i)).collect(Collectors.toList()));
    }

    /**
     * Saves copy-ratio-related posteriors to disk.
     *
     * TODO github/gatk-protected issue #855 -- write local copy ratio posteriors as well
     *
     * @param outputPath the output path
     */
    protected void writeCopyRatioPosteriors(final String outputPath) {
        logger.info("Saving copy ratio posteriors...");
        final List<List<HiddenStateSegmentRecord<STATE, Target>>> segments = getCopyRatioSegments();
        final String segmentsPath = new File(outputPath, CoverageModelGlobalConstants.COPY_RATIO_SEGMENTS_SUBDIR)
                .getAbsolutePath();
        createOutputPath(segmentsPath);

        sampleIndexStream().forEach(si -> {
            final File segmentsFile = new File(segmentsPath, processedSampleNameList.get(si) + ".seg");
            try (final HiddenStateSegmentRecordWriter<STATE, Target> segWriter =
                         new HiddenStateSegmentRecordWriter<>(segmentsFile)) {
                segWriter.writeAllRecords(segments.get(si));
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(segmentsFile, "Could not create copy ratio segments file");
            }
        });

        /* write the Viterbi copy number chains as a sample-target matrix */
        if (config.extendedPosteriorOutputEnabled()) {
            final File copyRatioViterbiFile = new File(outputPath, CoverageModelGlobalConstants.COPY_RATIO_VITERBI_FILENAME);
            final List<String> targetNames = processedReadCounts.targets().stream()
                    .map(Target::getName).collect(Collectors.toList());
            Nd4jIOUtils.writeNDArrayMatrixToTextFile(getViterbiAsNDArray(segments),
                    copyRatioViterbiFile, "VITERBI_COPY_RATIO_CHAIN", targetNames, processedSampleNameList);
        }

    }

    /**
     * Saves extended posteriors to disk
     *
     * @param outputPath the output path
     */
    protected void writeExtendedPosteriors(final String outputPath) {
        final List<String> sampleNames = processedReadCounts.columnNames();
        final List<String> targetNames = processedReadCounts.targets().stream()
                .map(Target::getName).collect(Collectors.toList());

        /* write total unexplained variance as a matrix */
        final File totalExplainedVarianceFile = new File(outputPath, CoverageModelGlobalConstants.TOTAL_UNEXPLAINED_VARIANCE_FILENAME);
        Nd4jIOUtils.writeNDArrayMatrixToTextFile(fetchTotalUnexplainedVariance(), totalExplainedVarianceFile,
                "TOTAL_UNEXPLAINED_VARIANCE", targetNames, sampleNames);

        if (biasCovariatesEnabled) {
            /* write total covariate bias per sample as a matrix */
            final File totalCovariateBiasFile = new File(outputPath, CoverageModelGlobalConstants.TOTAL_COVARIATE_BIAS_FILENAME);
            Nd4jIOUtils.writeNDArrayMatrixToTextFile(fetchTotalCovariateBiasPerSample(), totalCovariateBiasFile,
                    "TOTAL_COVARIATE_BIAS", targetNames, sampleNames);
        }
    }

    /**
     * Represents the verbosity level of output posteriors; see {@link #writePosteriors}
     */
    public enum PosteriorVerbosityLevel {
        /**
         * Basic posteriors (read depth, log likelihood, unexplained variance, bias covariates, ARD coefficients,
         * copy ratio segments, copy ratio max likelihood point estimates)
         */
        BASIC,

        /**
         * Extended posteriors (basic + total unexplained variance, total explained variance by bias covariates)
         */
        EXTENDED
    }
}
