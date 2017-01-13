package org.broadinstitute.hellbender.tools.coveragemodel;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.math3.analysis.UnivariateFunction;
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
import org.broadinstitute.hellbender.tools.coveragemodel.annotations.CachesRDD;
import org.broadinstitute.hellbender.tools.coveragemodel.annotations.EvaluatesRDD;
import org.broadinstitute.hellbender.tools.coveragemodel.annotations.UpdatesRDD;
import org.broadinstitute.hellbender.tools.coveragemodel.interfaces.CopyRatioExpectationsCalculator;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.FourierLinearOperator;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.FourierLinearOperatorNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.GeneralLinearOperator;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.IterativeLinearSolverNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.IterativeLinearSolverNDArray.ExitStatus;
import org.broadinstitute.hellbender.tools.coveragemodel.math.RobustBrentSolver;
import org.broadinstitute.hellbender.tools.coveragemodel.math.SynchronizedUnivariateSolver;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.GermlinePloidyAnnotatedTargetCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeData;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeDataCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.interfaces.AlleleMetadataProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.CallStringProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.ScalarProducer;
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
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelEMParams.CopyRatioHMMType.COPY_RATIO_HMM_LOCAL;

/**
 * This class represents the driver-node workspace for EM algorithm calculations of the coverage model
 *
 * @param <S> copy ratio (or copy number) state type
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class CoverageModelEMWorkspace<S extends AlleleMetadataProducer & CallStringProducer &
        ScalarProducer> {

    protected final Logger logger = LogManager.getLogger(CoverageModelEMWorkspace.class);

    protected final CoverageModelEMParams params;

    /**
     * The input read count collection after initial processing
     */
    protected final ReadCountCollection processedReadCounts;

    /**
     * Target list after initial processing
     */
    protected final List<Target> processedTargetList;

    /**
     * Same names list after initial processing
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
    protected final CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, S> copyRatioExpectationsCalculator;

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

    /* Spark-related members --- BEGIN */

    /**
     * The latest RDD of compute blocks
     */
    private JavaPairRDD<LinearSpaceBlock, CoverageModelEMComputeBlock> computeRDD;

    /**
     * The latest checkpointed RDD of compute blocks
     */
    private JavaPairRDD<LinearSpaceBlock, CoverageModelEMComputeBlock> prevCheckpointedComputeRDD;

    /**
     * A deque of cached RDDs of compute blocks
     */
    private Deque<JavaPairRDD<LinearSpaceBlock, CoverageModelEMComputeBlock>> prevCachedComputeRDDDeque = new LinkedList<>();

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
    protected List<LinearSpaceBlock> targetBlocks;

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
     * Norm_2 of bias covariates
     */
    protected final INDArray biasCovariatesNorm2;

    /**
     * Fourier factors of the CNV-avoiding regularizer
     */
    private final double[] fourierFactors;

    /**
     * Public constructor
     *
     * @param rawReadCounts an instance of {@link ReadCountCollection} containing raw read counts
     * @param germlinePloidyAnnotatedTargetCollection an instance of {@link GermlinePloidyAnnotatedTargetCollection}
     *                                                for obtaining target ploidies for different sex genotypes
     * @param sexGenotypeDataCollection an instance of {@link SexGenotypeDataCollection} for obtaining sample sex genotypes
     * @param params EM algorithm parameters
     * @param ctx the Spark context
     * @param copyRatioExpectationsCalculator an implementation of {@link CopyRatioExpectationsCalculator}
     */
    @UpdatesRDD @CachesRDD @EvaluatesRDD
    protected CoverageModelEMWorkspace(@Nonnull final ReadCountCollection rawReadCounts,
                                       @Nonnull final GermlinePloidyAnnotatedTargetCollection germlinePloidyAnnotatedTargetCollection,
                                       @Nonnull final SexGenotypeDataCollection sexGenotypeDataCollection,
                                       @Nonnull final CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, S> copyRatioExpectationsCalculator,
                                       @Nonnull final CoverageModelEMParams params,
                                       @Nullable final CoverageModelParameters model,
                                       @Nullable final JavaSparkContext ctx) {
        this.params = Utils.nonNull(params, "Coverage model EM-algorithm parameters must be non-null");
        this.copyRatioExpectationsCalculator = Utils.nonNull(copyRatioExpectationsCalculator, "Copy ratio posterior calculator" +
                " must be non-null");
        this.germlinePloidyAnnotatedTargetCollection = Utils.nonNull(germlinePloidyAnnotatedTargetCollection,
                "The germline ploidy-annotated target collection must be non-null");
        this.sexGenotypeDataCollection = Utils.nonNull(sexGenotypeDataCollection,
                "The sex genotype data collection must be non-null");
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
            final ReadCountCollection intermediateReadCounts = processReadCountCollection(targetSortedRawReadCounts,
                    params, logger);
            /* adapt model and read counts */
            final ImmutablePair<CoverageModelParameters, ReadCountCollection> modelReadCountsPair =
                    CoverageModelParameters.adaptModelToReadCountCollection(model, intermediateReadCounts, logger);
            processedModel = modelReadCountsPair.left;
            processedReadCounts = modelReadCountsPair.right;
            numLatents = processedModel.getNumLatents();
            if (params.getNumLatents() != processedModel.getNumLatents()) {
                logger.info("Changing number of latent variables to " + processedModel.getNumLatents() + " based" +
                        " on the provided model; requested value was: " + params.getNumLatents());
                params.setNumLatents(processedModel.getNumLatents());
            }
        } else {
            processedModel = null;
            processedReadCounts = processReadCountCollection(targetSortedRawReadCounts, params, logger);
            numLatents = params.getNumLatents();
        }

        numSamples = processedReadCounts.columnNames().size();
        numTargets = processedReadCounts.targets().size();
        processedTargetList = Collections.unmodifiableList(processedReadCounts.targets());
        processedSampleNameList = Collections.unmodifiableList(processedReadCounts.columnNames());
        processedSampleSexGenotypeData = Collections.unmodifiableList(processedSampleNameList.stream()
                .map(sexGenotypeDataCollection::getSampleSexGenotypeData)
                .collect(Collectors.toList()));
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
            this.numTargetBlocks = ParamUtils.inRange(params.getNumTargetSpaceParititions(), 1, numTargets,
                    "Number of target blocks must be between 1 and the size of target space.");
        }

        /* allocate memory for driver-node posterior expectations */
        sampleMeanLogReadDepths = Nd4j.zeros(numSamples, 1);
        sampleVarLogReadDepths = Nd4j.zeros(numSamples, 1);
        sampleBiasLatentPosteriorFirstMoments = Nd4j.zeros(numSamples, numLatents);
        sampleBiasLatentPosteriorSecondMoments = Nd4j.zeros(numSamples, numLatents, numLatents);
        sampleUnexplainedVariance = Nd4j.zeros(numSamples, 1);
        biasCovariatesNorm2 = Nd4j.zeros(1, numLatents);

        /* initialize the CNV-avoiding regularizer filter */
        if (params.fourierRegularizationEnabled()) {
            fourierFactors = FourierLinearOperator.getMidpassFilterFourierFactors(numTargets,
                    numTargets/params.getMaximumCNVLength(), numTargets/params.getMinimumCNVLength());
        } else {
            fourierFactors = null;
        }

        /* initialize target-space blocks */
        initializeTargetSpaceBlocks();

        /* create compute blocks */
        instantiateWorkers();

        /* push read counts to compute blocks */
        pushInitialDataToComputeBlocks();

        if (processedModel == null) {
            /* initialize model parameters and posterior expectations to default initial values */
            initializeModelParametersToDefaultValues();
        } else {
            initializeModelParametersFromGivenModel(processedModel);
        }
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
     * Process raw read counts and filter bad targets and/or samples as follows:
     *
     * <dl>
     *     <dt> Remove targets with uniformly low coverage </dt>
     * </dl>
     *
     * TODO github/gatk-protected issue #855 -- consider adding more filters
     *
     * - remove targets with very high and very low GC content (can be done externally)
     * - remove targets with lots of repeats (can be done externally)
     * - in the learning mode, remove a target if too many are masked across the samples (in that case, max likelihood
     *   parameter estimation is unreliable)
     *
     * @param rawReadCounts raw read counts
     * @return processed read counts
     */
    private static ReadCountCollection processReadCountCollection(@Nonnull final ReadCountCollection rawReadCounts,
                                                                 @Nonnull final CoverageModelEMParams params,
                                                                 @Nonnull final Logger logger) {
        ReadCountCollection processedReadCounts;
        processedReadCounts = ReadCountCollectionUtils.removeColumnsWithBadValues(rawReadCounts, logger);
        processedReadCounts = ReadCountCollectionUtils.removeTargetsWithUniformlyLowCoverage(processedReadCounts,
                params.getMinLearningReadCount(), logger);

        return processedReadCounts;
    }


    /**
     * Partitions the target space into {@link #numTargetBlocks} contiguous blocks
     */
    private void initializeTargetSpaceBlocks() {
        targetBlocks = CoverageModelSparkUtils.createLinearSpaceBlocks(numTargets, numTargetBlocks,
                CoverageModelGlobalConstants.DEFAULT_MIN_TARGET_BLOCK_SIZE);
        logger.debug("Target space blocks: " + targetBlocks.stream().map(LinearSpaceBlock::toString)
                .reduce((L, R) -> L + "\t" + R).orElse("None"));
    }

    /**
     * Determines whether a target with read count {@code readCount} and germline ploidy {@code germlinePloidy}
     * should be masked in the parameter estimation ("learning") step
     *
     * @param readCount raw read count
     * @param germlinePloidy germline ploidy
     * @return a boolean
     */
    private boolean isMaskedForLearning(final int readCount, final int germlinePloidy) {
        return readCount < params.getMinLearningReadCount() || germlinePloidy == 0;
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
        final List<ReadCountRecord> recs = processedReadCounts.records();

        logger.info("Pushing read count data to the worker(s)");
        final List<Tuple2<LinearSpaceBlock, CoverageModelEMComputeBlock.InitialDataBlock>> dataBlockList =
                new ArrayList<>();

        targetBlockStream().forEach(tb -> {
            /* take a contiguous [targets in the block x all samples] chunk from the read count collection
             * and ravel it in Fortran order */
            final int[] rawReadCountBlock = IntStream.range(tb.getBegIndex(), tb.getEndIndex())
                    /* 1-to-S flat map of each target to the read counts of all samples */
                    .mapToObj(ti -> recs.get(ti).getDoubleCounts())
                    .flatMapToDouble(Arrays::stream)
                    .mapToInt(d -> (int)FastMath.round(d)) /* round to integer values */
                    .toArray();

            /* fetch the germline ploidy within the same contiguous block of reads */
            final int[] germlinePloidyBlock = IntStream.range(tb.getBegIndex(), tb.getEndIndex())
                    /* map target index to actual targets */
                    .mapToObj(processedTargetList::get)
                    /* 1-to-S flat map of each target to the germline ploidy of all samples */
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
                    .forEach(idx -> germlinePloidyBlock[idx] =
                            CoverageModelGlobalConstants.READ_COUNT_ON_ZERO_PLOIDY_TARGETS);

            /* targets with low read count or zero ploidy are masked */
            final int[] maskBlock = IntStream.range(0, rawReadCountBlock.length)
                    .map(idx -> isMaskedForLearning(rawReadCountBlock[idx], germlinePloidyBlock[idx]) ? 0 : 1)
                    .toArray();

            /* TODO github/gatk-protected issue #855 -- in the future, this must be replaced with a sample-
             * and target-specific value calculated from the mapping quality distribution of each target */
            final double[] mappingErrorRateBlock = IntStream.range(0, rawReadCountBlock.length)
                    .mapToDouble(idx -> params.getMappingErrorRate())
                    .toArray();

            /* calculate copy ratio prior expectations */
            final int[] activeTargetIndices = IntStream.range(tb.getBegIndex(), tb.getEndIndex()).toArray();
            final List<Target> activeTargets = Arrays.stream(activeTargetIndices)
                    .mapToObj(processedTargetList::get)
                    .collect(Collectors.toList());
            final List<CopyRatioExpectations> copyRatioPriorExpectationsList = sampleIndexStream()
                    .mapToObj(si -> copyRatioExpectationsCalculator.getCopyRatioPriorExpectations(
                            CopyRatioCallingMetadata.builder()
                                    .setSampleIndex(si)
                                    .setSampleName(processedSampleNameList.get(si))
                                    .setSampleSexGenotypeData(processedSampleSexGenotypeData.get(si))
                                    .build(), activeTargets))
                    .collect(Collectors.toList());

            final double[] logCopyRatioPriorMeansBlock = IntStream.range(0, tb.getNumTargets())
                    .mapToObj(rti -> copyRatioPriorExpectationsList.stream()
                            .mapToDouble(cre -> cre.getLogCopyRatioMeans()[rti])
                            .toArray())
                    .flatMapToDouble(Arrays::stream)
                    .toArray();

            final double[] logCopyRatioPriorVariancesBlock = IntStream.range(0, tb.getNumTargets())
                    .mapToObj(rti -> copyRatioPriorExpectationsList.stream()
                            .mapToDouble(cre -> cre.getLogCopyRatioVariances()[rti])
                            .toArray())
                    .flatMapToDouble(Arrays::stream)
                    .toArray();

            /* we do not need to take care of log copy ratio means and variances on masked targets here.
             * potential NaNs will be rectified in the compute blocks by calling the method
             * {@link CoverageModelEMComputeBlock#cloneWithInitializedData} */

            /* add the block to list */
            dataBlockList.add(new Tuple2<>(tb, new CoverageModelEMComputeBlock.InitialDataBlock(rawReadCountBlock,
                    maskBlock, logCopyRatioPriorMeansBlock, logCopyRatioPriorVariancesBlock, mappingErrorRateBlock)));
        });

        /* push to compute blocks */
        joinWithWorkersAndMap(dataBlockList, p -> p._1.cloneWithInitializedData(p._2));
    }

    /**
     * Initialize model parameters to default values
     */
    @UpdatesRDD
    private void initializeModelParametersToDefaultValues() {
        /* make local references for lambda captures */
        final CoverageModelEMParams params = this.params;
        final int numSamples = this.numSamples;
        final int numLatents = params.getNumLatents();

        /* set $W_{tm}$ to a zero-padded D X D scaled identity matrix
         * set $m_t$ to zero
         * set $\Psi_t$ to zero
         * set $z_sl$ and $zz_sll$ to zero */
        mapWorkers(cb -> {
            /* generate the portion of the truncated identity that belongs to this target-space block */
            final LinearSpaceBlock tb = cb.getTargetSpaceBlock();
            final INDArray newBiasCovariates = Nd4j.zeros(tb.getNumTargets(), numLatents);
            if (tb.getBegIndex() < numLatents) {
                IntStream.range(tb.getBegIndex(), FastMath.min(tb.getEndIndex(), numLatents)).forEach(ti ->
                        newBiasCovariates.getRow(ti).assign(Nd4j.zeros(1, numLatents).putScalar(0, ti,
                                CoverageModelGlobalConstants.INITIAL_BIAS_COVARIATES_SCALAR)));
            }
            return cb.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl,
                    newBiasCovariates)
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t,
                            Nd4j.zeros(1, cb.getTargetSpaceBlock().getNumTargets()))
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t,
                            Nd4j.ones(1, cb.getTargetSpaceBlock().getNumTargets())
                            .mul(CoverageModelGlobalConstants.INITIAL_TARGET_UNEXPLAINED_VARIANCE))
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.z_sl,
                            Nd4j.zeros(numSamples, numLatents))
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.zz_sll,
                            Nd4j.zeros(numSamples, numLatents, numLatents))
                    .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.gamma_s,
                            Nd4j.zeros(numSamples, 1));
        });
        if (params.fourierRegularizationEnabled()) {
            mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock
                    .CoverageModelICGCacheTag.M_STEP_W_REG));
            updateFilteredBiasCovariates();
        }
    }

    /**
     * Initialize the compute blocks with given model parameters
     * @param model coverage model parameters
     */
    @UpdatesRDD
    private void initializeModelParametersFromGivenModel(@Nonnull final CoverageModelParameters model) {
        /* make local references for lambda capture */
        final CoverageModelEMParams params = this.params;
        final int numSamples = this.numSamples;
        final int numLatents = params.getNumLatents();

        /* sample-specific unexplained variances can not be lower than the following value */
        final double gammaLowerBound = 0.0;

        if (sparkContextIsAvailable) { /* broadcast the model */
            final Broadcast<CoverageModelParameters> broadcastedModel = ctx.broadcast(model);
            mapWorkers(cb -> {
                final LinearSpaceBlock tb = cb.getTargetSpaceBlock();
                return cb.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl,
                        broadcastedModel.getValue().getBiasCovariatesOnTargetBlock(tb))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t,
                                broadcastedModel.getValue().getTargetMeanBiasOnTargetBlock(tb))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t,
                                broadcastedModel.getValue().getTargetUnexplainedVarianceOnTargetBlock(tb))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.z_sl,
                                Nd4j.zeros(numSamples, numLatents))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.zz_sll,
                                Nd4j.zeros(numSamples, numLatents, numLatents))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.gamma_s,
                                Nd4j.ones(numSamples, 1).muli(gammaLowerBound));
            });
        } else { /* pass the model by reference */
            mapWorkers(cb -> {
                final LinearSpaceBlock tb = cb.getTargetSpaceBlock();
                return cb.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl,
                        model.getBiasCovariatesOnTargetBlock(tb))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t,
                                model.getTargetMeanBiasOnTargetBlock(tb))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t,
                                model.getTargetUnexplainedVarianceOnTargetBlock(tb))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.z_sl,
                                Nd4j.zeros(numSamples, numLatents))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.zz_sll,
                                Nd4j.zeros(numSamples, numLatents, numLatents))
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.gamma_s,
                                Nd4j.ones(numSamples, 1).muli(gammaLowerBound));
            });
        }

        /* update the local copy of gamma as well */
        sampleUnexplainedVariance.assign(Nd4j.ones(numSamples, 1).muli(gammaLowerBound));

        if (params.fourierRegularizationEnabled()) {
            mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.M_STEP_W_REG));
            updateFilteredBiasCovariates();
        }
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
        IntStream.range(0, numLatents).parallel()
                .forEach(li -> {
                    final INDArrayIndex[] slice = {NDArrayIndex.all(), NDArrayIndex.point(li)};
                    filteredBiasCovariates.get(slice).assign(
                            regularizerFourierLinearOperator.operate(biasCovariates.get(slice)));
                });

        /* sent the new W to workers */
        switch (params.getBiasCovariatesComputeNodeCommunicationPolicy()) {
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
        final double fact = params.getFourierRegularizationStrength() * psiAverage;
        return new FourierLinearOperatorNDArray(numTargets,
                Arrays.stream(fourierFactors).map(f -> f * fact).toArray(), params.zeroPadFFT());
    }


    /* E-step methods */

    /**
     * Updates the first (E[z]) and second (E[z z^T]) posterior moments of the log bias continuous
     * latent variables (z)
     *
     * @implNote the operations done on the driver node have low complexity only if D, the dimension of the latent
     * space, is small:
     *
     *     (a) G_s = (I + [contribGMatrix])^{-1} for each sample \sim O(S x D^3)
     *     (b) E[z_s] = G_s [contribZ_s] for each sample \sim O(S x D^3)
     *     (c) E[z_s z_s^T] = G_s + E[z_s] E[z_s^T] for each sample \sim O(S x D^2)
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateBiasLatentPosteriorExpectations() {
        /* calculate and cache the required quantities on compute blocks */
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_Z));
        cacheWorkers("after E-step for bias initialization");

        final ImmutableTriple<INDArray, INDArray, INDArray> biasLatentPosteriorExpectationsData =
                fetchBiasLatentPosteriorExpectationsDataFromWorkers();
        final INDArray contribGMatrix = biasLatentPosteriorExpectationsData.left;
        final INDArray contribZ = biasLatentPosteriorExpectationsData.middle;
        final INDArray sharedPart = biasLatentPosteriorExpectationsData.right;

        /* calculate G_{s\mu\nu} = (sharedPart + W^T M \Psi^{-1} W)^{-1} by doing sample-wise matrix inversion */
        final INDArray sampleGTensor = Nd4j.create(numSamples, numLatents, numLatents);
        sampleIndexStream().forEach(si -> sampleGTensor.get(NDArrayIndex.point(si)).assign(
                CoverageModelEMWorkspaceMathUtils.minv(sharedPart.add(contribGMatrix.get(NDArrayIndex.point(si))))));

        final INDArray newSampleBiasLatentPosteriorFirstMoments = Nd4j.create(numSamples, numLatents);
        final INDArray newSampleBiasLatentPosteriorSecondMoments = Nd4j.create(numSamples, numLatents, numLatents);

        sampleIndexStream().forEach(si -> {
            final INDArray sampleGMatrix = sampleGTensor.get(NDArrayIndex.point(si), NDArrayIndex.all(),
                    NDArrayIndex.all());
            /* E[z_s] = G_s W^T M_{st} \Psi_{st}^{-1} (m_{st} - m_t) */
            newSampleBiasLatentPosteriorFirstMoments.get(NDArrayIndex.point(si), NDArrayIndex.all())
                    .assign(sampleGMatrix.mmul(contribZ.get(NDArrayIndex.all(), NDArrayIndex.point(si))).transpose());
            /* E[z_s z_s^T] = G_s + E[z_s] E[z_s^T] */
            newSampleBiasLatentPosteriorSecondMoments.get(NDArrayIndex.point(si), NDArrayIndex.all(), NDArrayIndex.all())
                    .assign(sampleGMatrix.add(
                            newSampleBiasLatentPosteriorFirstMoments.get(NDArrayIndex.point(si), NDArrayIndex.all()).transpose()
                                    .mmul(newSampleBiasLatentPosteriorFirstMoments.get(NDArrayIndex.point(si), NDArrayIndex.all()))));
        });

        /* admix with old posteriors */
        final INDArray newSampleBiasLatentPosteriorFirstMomentsAdmixed = newSampleBiasLatentPosteriorFirstMoments
                .mul(params.getMeanFieldAdmixingRatio()).addi(sampleBiasLatentPosteriorFirstMoments
                        .mul(1.0 - params.getMeanFieldAdmixingRatio()));
        final INDArray newSampleBiasLatentPosteriorSecondMomentsAdmixed = newSampleBiasLatentPosteriorSecondMoments
                .mul(params.getMeanFieldAdmixingRatio()).addi(sampleBiasLatentPosteriorSecondMoments
                        .mul(1.0 - params.getMeanFieldAdmixingRatio()));

        /* calculate the error from the change in E[z_s] */
        final double errorNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                newSampleBiasLatentPosteriorFirstMomentsAdmixed.sub(sampleBiasLatentPosteriorFirstMoments));

        /* update driver-node copies */
        sampleBiasLatentPosteriorFirstMoments.assign(newSampleBiasLatentPosteriorFirstMomentsAdmixed);
        sampleBiasLatentPosteriorSecondMoments.assign(newSampleBiasLatentPosteriorSecondMomentsAdmixed);

        /* broadcast the new bias latent posteriors */
        pushToWorkers(ImmutablePair.of(newSampleBiasLatentPosteriorFirstMomentsAdmixed, newSampleBiasLatentPosteriorSecondMomentsAdmixed),
                (dat, cb) -> cb
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.z_sl, dat.left)
                        .cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.zz_sll, dat.right));

        return SubroutineSignal.builder().put("error_norm", errorNormInfinity).build();
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
        if (params.fourierRegularizationEnabled()) {
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
    public SubroutineSignal updateReadDepthPosteriorExpectations() {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_D));
        cacheWorkers("after E-step for read depth initialization");
        /* map each compute block to their respective read depth estimation data (a triple of rank-1 INDArray's each
         * with S elements) and reduce by pairwise addition */
        final ImmutablePair<INDArray, INDArray> factors = mapWorkersAndReduce(
                CoverageModelEMComputeBlock::getReadDepthLatentPosteriorData,
                (p1, p2) -> ImmutablePair.of(p1.left.add(p2.left), p1.right.add(p2.right)));

        /* put together */
        final INDArray numerator = factors.left;
        final INDArray denominator = factors.right;

        final INDArray newSampleMeanLogReadDepths = numerator.div(denominator);
        final INDArray newSampleVarLogReadDepths = Nd4j.ones(denominator.shape()).div(denominator);

        /* admix */
        final INDArray newSampleMeanLogReadDepthsAdmixed = newSampleMeanLogReadDepths
                .mul(params.getMeanFieldAdmixingRatio())
                .addi(sampleMeanLogReadDepths.mul(1.0 - params.getMeanFieldAdmixingRatio()));
        final INDArray newSampleVarLogReadDepthsAdmixed = newSampleVarLogReadDepths
                .mul(params.getMeanFieldAdmixingRatio())
                .addi(sampleVarLogReadDepths.mul(1.0 - params.getMeanFieldAdmixingRatio()));

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

        return SubroutineSignal.builder().put("error_norm", errorNormInfinity).build();
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

        final double gammaLowerLimit = 0.0;
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
                    cb.calculateGammaObjectiveFunctionMultiSample(sampleIndices, gammaValues), INDArray::add);
            /* create output map */
            final Map<Integer, Double> output = new HashMap<>();
            IntStream.range(0, sampleIndices.length)
                    .forEach(evalIdx -> output.put(sampleIndices[evalIdx], eval.getDouble(evalIdx)));
            return output;
        };

        /* instantiate a synchronized multi-sample root finder and add jobs */
        final SynchronizedUnivariateSolver syncSolver = new SynchronizedUnivariateSolver(objFunc,
                numSamples, RobustBrentSolver.MeritPolicy.LARGEST_ROOT, params.getGammaSolverNumBisections(),
                params.getGammaSolverRefinementDepth());
        IntStream.range(0, numSamples)
                .forEach(si -> {
                    final double x0 = 0.5 * (gammaLowerLimit + params.getGammaUpperLimit());
                    syncSolver.add(si, gammaLowerLimit, params.getGammaUpperLimit(), x0,
                            params.getGammaAbsoluteTolerance(), params.getGammaRelativeTolerance(),
                            params.getGammaMaximumIterations());
                });

        /* solve and collect statistics */
        final INDArray newSampleUnexplainedVariance = Nd4j.create(numSamples, 1);
        final List<Integer> numberOfEvaluations = new ArrayList<>(numSamples);
        try {
            final Map<Integer, SynchronizedUnivariateSolver.UnivariateSolverSummary> newGammaMap = syncSolver.solve();
            newGammaMap.entrySet().forEach(entry -> {
                final int sampleIndex = entry.getKey();
                final SynchronizedUnivariateSolver.UnivariateSolverSummary summary = entry.getValue();
                double val =  gammaLowerLimit;
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
                .mul(params.getMeanFieldAdmixingRatio())
                .addi(sampleUnexplainedVariance.mul(1 - params.getMeanFieldAdmixingRatio()));

        /* calculate the error */
        final double errorNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                newSampleUnexplainedVarianceAdmixed.sub(sampleUnexplainedVariance));

        /* update local copy */
        sampleUnexplainedVariance.assign(newSampleUnexplainedVarianceAdmixed);

        /* push to workers */
        pushToWorkers(newSampleUnexplainedVarianceAdmixed, (arr, cb) -> cb.cloneWithUpdatedPrimitive(
                CoverageModelEMComputeBlock.CoverageModelICGCacheNode.gamma_s,
                newSampleUnexplainedVarianceAdmixed));

        return SubroutineSignal.builder()
                .put("error_norm", errorNormInfinity)
                .put("iterations", (int)(numberOfEvaluations.stream().mapToDouble(d -> d).sum() / numSamples))
                .build();
    }

    /**
     * E-step update of copy ratio posteriors
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    public SubroutineSignal updateCopyRatioPosteriorExpectations() {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_C));
        cacheWorkers("after E-step for copy ratio initialization");

        /* calculate posteriors */
        logger.debug("Copy ratio HMM type: " + params.getCopyRatioHMMType().name());
        long startTime = System.nanoTime();
        final SubroutineSignal sig;
        if (params.getCopyRatioHMMType().equals(COPY_RATIO_HMM_LOCAL) || !sparkContextIsAvailable) {
            /* local mode */
            sig = updateCopyRatioPosteriorExpectationsLocal();
        } else {
            /* spark mode */
            sig = updateCopyRatioPosteriorExpectationsSpark();
        }
        long endTime = System.nanoTime();
        logger.debug("Copy ratio posteriors calculation time: " + (double)(endTime - startTime)/1000000 + " ms");
        return sig;
    }

    /**
     * Fetches forward-backward and Viterbi algorithm results on all samples
     *
     * @return a list of {@link CopyRatioHiddenMarkovModelResults}
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    protected List<CopyRatioHiddenMarkovModelResults<CoverageModelCopyRatioEmissionData, S>> getCopyRatioHiddenMarkovModelResults() {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.E_STEP_C));
        cacheWorkers("after E-step for copy ratio HMM result generation");
        final List<CopyRatioHiddenMarkovModelResults<CoverageModelCopyRatioEmissionData, S>> result;
        /* calculate posteriors */
        logger.debug("Copy ratio HMM type: " + params.getCopyRatioHMMType().name());
        long startTime = System.nanoTime();
        if (params.getCopyRatioHMMType().equals(COPY_RATIO_HMM_LOCAL) || !sparkContextIsAvailable) {
            /* local mode */
            result = getCopyRatioHiddenMarkovModelResultsLocal();
        } else {
            /* spark mode */
            result = getCopyRatioHiddenMarkovModelResultsSpark();
        }
        long endTime = System.nanoTime();
        logger.debug("Copy ratio HMM result generation time: " + (double)(endTime - startTime)/1000000 + " ms");
        return result;
    }

    /**
     * Local implementation of the E-step update of copy ratio posteriors
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    public SubroutineSignal updateCopyRatioPosteriorExpectationsLocal() {
        /* step 1. fetch copy ratio emission data */
        final List<List<CoverageModelCopyRatioEmissionData>> copyRatioEmissionData = fetchCopyRatioEmissionDataLocal();

        /* step 2. run the forward-backward algorithm and calculate copy ratio posteriors */
        final INDArray sampleReadDepths = Transforms.exp(sampleMeanLogReadDepths, true);
        final List<CopyRatioExpectations> copyRatioPosteriorResults = sampleIndexStream()
                .parallel()
                .mapToObj(si -> copyRatioExpectationsCalculator.getCopyRatioPosteriorExpectations(
                        CopyRatioCallingMetadata.builder()
                                .setSampleIndex(si)
                                .setSampleName(processedSampleNameList.get(si))
                                .setSampleSexGenotypeData(processedSampleSexGenotypeData.get(si))
                                .setSampleCoverageDepth(sampleReadDepths.getDouble(si))
                                .setEmissionCalculationStrategy(CopyRatioCallingMetadata.EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                                .build(),
                        processedTargetList,
                        copyRatioEmissionData.get(si)))
                .collect(Collectors.toList());

        /* sent the results back to workers */
        final ImmutablePair<INDArray, INDArray> copyRatioPosteriorDataPair =
                convertCopyRatioLatentPosteriorExpectationsToNDArray(copyRatioPosteriorResults);
        final INDArray log_c_st = copyRatioPosteriorDataPair.left;
        final INDArray var_log_c_st = copyRatioPosteriorDataPair.right;

        final double meanFieldAdmixingRatio = params.getMeanFieldAdmixingRatio();

        /* partition the pair of (log_c_st, var_log_c_st), sent the result to workers via broadcast-hash-map */
        pushToWorkers(mapINDArrayPairToBlocks(log_c_st.transpose(), var_log_c_st.transpose()),
                (p, cb) -> cb.cloneWithUpdatedCopyRatioPosteriors(
                        p.get(cb.getTargetSpaceBlock()).left.transpose(),
                        p.get(cb.getTargetSpaceBlock()).right.transpose(),
                        meanFieldAdmixingRatio));
        cacheWorkers("after E-step update of copy ratio posteriors");

        /* collect subroutine signals */
        final List<SubroutineSignal> sigs = mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal);

        final double errorNormInfinity = Collections.max(sigs.stream()
                .map(sig -> sig.getDouble("error_norm"))
                .collect(Collectors.toList()));

        return SubroutineSignal.builder()
                .put("error_norm", errorNormInfinity)
                .build();
    }

    /**
     * Queries copy ratio emission data from compute blocks
     *
     * @return a double list of {@link CoverageModelCopyRatioEmissionData}
     */
    private List<List<CoverageModelCopyRatioEmissionData>> fetchCopyRatioEmissionDataLocal() {
        /* fetch data from workers */
        final List<ImmutablePair<LinearSpaceBlock, List<List<CoverageModelCopyRatioEmissionData>>>> collectedCopyRatioData =
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

        sampleIndexStream().forEach(si -> {
            final CopyRatioExpectations res = copyRatioExpectationsList.get(si);
            log_c_st.getRow(si).assign(Nd4j.create(res.getLogCopyRatioMeans(), new int[] {1, numTargets}));
            var_log_c_st.getRow(si).assign(Nd4j.create(res.getLogCopyRatioVariances(), new int[] {1, numTargets}));
        });
        return ImmutablePair.of(log_c_st, var_log_c_st);
    }

    private List<CopyRatioHiddenMarkovModelResults<CoverageModelCopyRatioEmissionData, S>> getCopyRatioHiddenMarkovModelResultsLocal() {
        final List<List<CoverageModelCopyRatioEmissionData>> copyRatioEmissionData = fetchCopyRatioEmissionDataLocal();
        final INDArray sampleReadDepths = Transforms.exp(sampleMeanLogReadDepths, true);
        return sampleIndexStream()
                .mapToObj(si -> copyRatioExpectationsCalculator.getCopyRatioHiddenMarkovModelResults(
                        CopyRatioCallingMetadata.builder()
                                .setSampleIndex(si)
                                .setSampleName(processedSampleNameList.get(si))
                                .setSampleSexGenotypeData(processedSampleSexGenotypeData.get(si))
                                .setSampleCoverageDepth(sampleReadDepths.getDouble(si))
                                .setEmissionCalculationStrategy(CopyRatioCallingMetadata.EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                                .build(),
                        processedTargetList,
                        copyRatioEmissionData.get(si)))
                .collect(Collectors.toList());
    }

    /**
     * The Spark implementation of the E-step update of copy ratio posteriors
     *
     * @return a {@link SubroutineSignal} containing the update size
     */
    @EvaluatesRDD @UpdatesRDD @CachesRDD
    private SubroutineSignal updateCopyRatioPosteriorExpectationsSpark() {
        /* local final member variables for lambda capture */
        final List<LinearSpaceBlock> targetBlocks = new ArrayList<>();
        targetBlocks.addAll(this.targetBlocks);
        final List<Target> targetList = new ArrayList<>();
        targetList.addAll(processedTargetList);
        final List<String> sampleNameList = new ArrayList<>();
        sampleNameList.addAll(processedSampleNameList);
        final List<SexGenotypeData> sampleSexGenotypeData = new ArrayList<>();
        sampleSexGenotypeData.addAll(processedSampleSexGenotypeData);
        final int numTargetBlocks = targetBlocks.size();
        final CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, S> calculator =
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
                                        .setSampleIndex(si)
                                        .setSampleName(sampleNameList.get(si))
                                        .setSampleSexGenotypeData(sampleSexGenotypeData.get(si))
                                        .setSampleCoverageDepth(sampleReadDepths.getDouble(si))
                                        .setEmissionCalculationStrategy(CopyRatioCallingMetadata.EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                                        .build();
                                newPartitionData.add(new Tuple2<>(prevDatum._1,
                                        calculator.getCopyRatioPosteriorExpectations(copyRatioCallingMetadata,
                                                targetList, prevDatum._2)));
                            }
                            return newPartitionData.iterator();
                        }, true);

        /* repartition in target space */
        final JavaPairRDD<LinearSpaceBlock, ImmutablePair<INDArray, INDArray>>
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

        /* step 3. merge with computeRDD and update */
        final double admixingRatio = params.getMeanFieldAdmixingRatio();
        computeRDD = computeRDD.join(blockifiedCopyRatioPosteriorResultsPairRDD)
                .mapValues(t -> t._1.cloneWithUpdatedCopyRatioPosteriors(t._2.left, t._2.right, admixingRatio));
        cacheWorkers("after E-step for copy ratio update");

        /* collect subroutine signals */
        final List<SubroutineSignal> sigs = mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal);

        final double errorNormInfinity = Collections.max(sigs.stream()
                .map(sig -> sig.getDouble("error_norm"))
                .collect(Collectors.toList()));

        return SubroutineSignal.builder()
                .put("error_norm", errorNormInfinity)
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
                    final LinearSpaceBlock tb = tuple._1;
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
                        /* remove the LinearSpaceBlock keys from the sorted emissionBlocksList */
                        .map(p -> p._2)
                        /* flatten */
                        .flatMap(List::stream)
                        /* collect as a single list */
                        .collect(Collectors.toList()));
    }

    /**
     * Stack a list of (1 x T INDArray, 1 x T INDArray) pairs along the 0th axis of each INDArray
     * and returns a (S x T INDArray, S x T INDArray) pair
     *
     * @param perSampleData a list of (1 x T INDArray, 1 x T INDArray) pairs
     * @return a (S x T INDArray, S x T INDArray) pair
     */
    private static ImmutablePair<INDArray, INDArray> stackCopyRatioPosteriorDataForAllSamples(
            final List<ImmutablePair<INDArray, INDArray>> perSampleData) {
        return ImmutablePair.of(Nd4j.vstack(perSampleData.stream().map(p -> p.left).collect(Collectors.toList())),
                Nd4j.vstack(perSampleData.stream().map(p -> p.right).collect(Collectors.toList())));
    }

    /**
     * Fetch forward-backward and Viterbi results from compute blocks (Spark implementation)
     *
     * @return a list of {@link CopyRatioHiddenMarkovModelResults}
     */
    private List<CopyRatioHiddenMarkovModelResults<CoverageModelCopyRatioEmissionData, S>> getCopyRatioHiddenMarkovModelResultsSpark() {
        /* local final member variables for lambda capture */
        final List<Target> targetList = new ArrayList<>();
        targetList.addAll(processedTargetList);
        final List<SexGenotypeData> sampleSexGenotypeData = new ArrayList<>();
        sampleSexGenotypeData.addAll(processedSampleSexGenotypeData);
        final List<String> sampleNameList = new ArrayList<>();
        sampleNameList.addAll(processedSampleNameList);
        final INDArray sampleReadDepths = Transforms.exp(sampleMeanLogReadDepths, true);
        final CopyRatioExpectationsCalculator<CoverageModelCopyRatioEmissionData, S> calculator =
                this.copyRatioExpectationsCalculator;

        return fetchCopyRatioEmissionDataSpark()
                /* let the workers run fb and Viterbi */
                .mapPartitionsToPair(it -> {
                    final List<Tuple2<Integer, CopyRatioHiddenMarkovModelResults<
                            CoverageModelCopyRatioEmissionData, S>>> newPartitionData = new ArrayList<>();
                    while (it.hasNext()) {
                        final Tuple2<Integer, List<CoverageModelCopyRatioEmissionData>> prevDatum = it.next();
                        final int sampleIndex = prevDatum._1;
                        final CopyRatioCallingMetadata copyRatioCallingMetadata = CopyRatioCallingMetadata.builder()
                                .setSampleIndex(sampleIndex)
                                .setSampleName(sampleNameList.get(sampleIndex))
                                .setSampleSexGenotypeData(sampleSexGenotypeData.get(sampleIndex))
                                .setSampleCoverageDepth(sampleReadDepths.getDouble(sampleIndex))
                                .setEmissionCalculationStrategy(CopyRatioCallingMetadata.EmissionCalculationStrategy.HYBRID_POISSON_GAUSSIAN)
                                .build();
                        newPartitionData.add(new Tuple2<>(sampleIndex,
                                calculator.getCopyRatioHiddenMarkovModelResults(copyRatioCallingMetadata,
                                        targetList, prevDatum._2)));
                    }
                    return newPartitionData.iterator();
                }, true)
                /* collect the data to driver node a list of [sample index, HMM results] */
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
                mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal)
                        .stream().map(sig -> sig.getDouble("error_norm")).collect(Collectors.toList()));

        return SubroutineSignal.builder().put("error_norm", errorNormInfinity).build();
    }

    /**
     * M-step for target unexplained variance ($\Psi_t$)
     *
     * @return a {@link SubroutineSignal} object containing information about the solution
     */
    @UpdatesRDD @CachesRDD
    public SubroutineSignal updateTargetUnexplainedVariance() {
        final int psiMaxIterations = params.getPsiMaxIterations();
        final double psiAbsoluteTolerance = params.getPsiAbsoluteTolerance();
        final double psiRelativeTolerance = params.getPsiRelativeTolerance();
        final double psiUpperLimit = params.getPsiUpperLimit();
        final int psiSolverNumBisections = params.getPsiSolverNumBisections();
        final int psiSolverRefinementDepth = params.getPsiSolverRefinementDepth();

        logger.debug("Psi solver type: " + params.getPsiUpdateMode().name());
        switch (params.getPsiUpdateMode()) {
            case PSI_TARGET_RESOLVED: /* done on the compute blocks */
                mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.M_STEP_PSI)
                        .cloneWithUpdatedTargetUnexplainedVarianceTargetResolved(psiMaxIterations, psiUpperLimit, psiAbsoluteTolerance,
                                psiRelativeTolerance, psiSolverNumBisections, psiSolverRefinementDepth));
                break;

            case PSI_ISOTROPIC: /* done on the driver node */
                return updateTargetUnexplainedVarianceIsotropic();

            default:
                throw new RuntimeException("Illegal Psi solver type.");
        }

        cacheWorkers("after M-step for target unexplained variance");

        /* accumulate error from all workers */
        final List<SubroutineSignal> signalList = mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal);
        final double errorNormInfinity = Collections.max(signalList.stream().map(sig -> sig.getDouble("error_norm"))
                .collect(Collectors.toList()));
        final int maxIterations = Collections.max(signalList.stream().map(sig -> sig.getInteger("iterations"))
                .collect(Collectors.toList()));
        final int minIterations = Collections.min(signalList.stream().map(sig -> sig.getInteger("iterations"))
                .collect(Collectors.toList()));
        return SubroutineSignal.builder()
                .put("error_norm", errorNormInfinity)
                .put("min_iterations", minIterations)
                .put("max_iterations", maxIterations)
                .put("iterations", maxIterations) /* for uniformity */
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

        final double oldIsotropicPsi = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t, 1)
                .meanNumber().doubleValue();

        final double psiLowerBound = 0.0;
        final UnivariateFunction objFunc = psi -> mapWorkersAndReduce(cb ->
                cb.calculateSampleTargetSummedPsiObjectiveFunction(psi), (a, b) -> a + b);
        final UnivariateFunction meritFunc = psi -> mapWorkersAndReduce(cb ->
                cb.calculateSampleTargetSummedPsiMeritFunction(psi), (a, b) -> a + b);

        final RobustBrentSolver solver = new RobustBrentSolver(params.getPsiRelativeTolerance(),
                params.getPsiAbsoluteTolerance(), CoverageModelGlobalConstants.DEFAULT_FUNCTION_EVALUATION_ACCURACY);
        double newIsotropicPsi;
        try {
            newIsotropicPsi = solver.solve(params.getPsiMaxIterations(), objFunc, meritFunc, null,
                    psiLowerBound, params.getPsiUpperLimit(), params.getPsiSolverNumBisections(),
                    params.getPsiSolverRefinementDepth());
        } catch (NoBracketingException e) {
            logger.warn("Root of M-step optimality equation for isotropic unexplained variance could be bracketed");
            newIsotropicPsi = oldIsotropicPsi;
        } catch (TooManyEvaluationsException e) {
            logger.warn("Too many evaluations -- increase the number of root-finding iterations for the M-step update" +
                    " of unexplained variance");
            newIsotropicPsi = oldIsotropicPsi;
        }

        /* update the compute block(s) */
        final double errNormInfinity = FastMath.abs(newIsotropicPsi - oldIsotropicPsi);
        final int maxIterations = solver.getEvaluations();
        final double finalizedNewIsotropicPsi = newIsotropicPsi;
        mapWorkers(cb -> cb.cloneWithUpdatedPrimitive(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Psi_t,
                Nd4j.ones(1, cb.getTargetSpaceBlock().getNumTargets()).muli(finalizedNewIsotropicPsi)));
        return SubroutineSignal.builder()
                .put("error_norm", errNormInfinity)
                .put("iterations", maxIterations).build();
    }

    /**
     * M-step update of bias covariates
     *
     * @return a {@link SubroutineSignal} object containing "error_norm"
     */
    @UpdatesRDD @CachesRDD
    public SubroutineSignal updateBiasCovariates() {
        /* perform the M-step update */
        final SubroutineSignal sig;
        if (!params.fourierRegularizationEnabled()) {
            sig = updateBiasCovariatesUnregularized();
        } else {
            sig = updateBiasCovariatesRegularized();
        }
        return sig;
    }

    /**
     * M-step update of bias covariates w/o regularization
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    @UpdatesRDD @CachesRDD
    private SubroutineSignal updateBiasCovariatesUnregularized() {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.M_STEP_W_UNREG)
                .cloneWithUpdatedBiasCovariatesUnregularized());
        cacheWorkers("after M-step update of bias covariates w/o regularization");

        /* orthogonalize if required */
        if (params.isOrthogonalizeAndSortBiasCovariatesEnabled()) {
            orthogonalizeAndSortBiasCovariates();
            cacheWorkers("after orthogonalization of bias covariates");
        }

        /* accumulate error from all nodes */
        final double errorNormInfinity = Collections.max(
                mapWorkersAndCollect(CoverageModelEMComputeBlock::getLatestMStepSignal)
                        .stream().map(sig -> sig.getDouble("error_norm")).collect(Collectors.toList()));
        return SubroutineSignal.builder().put("error_norm", errorNormInfinity).build();
    }

    /**
     * M-step update of bias covariates w/ regularization (local implementation)
     *
     * @return a {@link SubroutineSignal} containing the update size (key: "error_norm")
     */
    @UpdatesRDD @EvaluatesRDD @CachesRDD
    private SubroutineSignal updateBiasCovariatesRegularized() {
        mapWorkers(cb -> cb.cloneWithUpdatedCachesByTag(CoverageModelEMComputeBlock.CoverageModelICGCacheTag.M_STEP_W_REG));
        cacheWorkers("after M-step update of bias covariates w/ regularization");

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
                params.getWAbsoluteTolerance(), params.getWRelativeTolerance(), params.getWMaxIterations(),
                x -> x.normmaxNumber().doubleValue(), /* norm */
                (x, y) -> x.mul(y).sumNumber().doubleValue(), /* inner product */
                true);

        /* solve */
        long startTime = System.nanoTime();
        final SubroutineSignal sig = iterSolver.cg(W_tl_old);
        linop.cleanupAfter();
        precond.cleanupAfter();
        long endTime = System.nanoTime();
        logger.debug("CG execution time for solving the regularized M-step update equation for bias covariates" +
                (double)(endTime - startTime)/1000000 + " ms");

        /* check the exit status of the solver and push the new W to workers */
        final ExitStatus exitStatus = (ExitStatus)sig.getObject("status");
        if (exitStatus == ExitStatus.FAIL_MAX_ITERS) {
            logger.warn("CG iterations for M-step update of bias covariates did not converge. Increase maximum iterations" +
                    " and/or decrease absolute/relative error tolerances");
        }
        final int iters = sig.getInteger("iterations");
        final INDArray W_tl_new = sig.getINDArray("x");

        switch (params.getBiasCovariatesComputeNodeCommunicationPolicy()) {
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

        /* orthogonalize if required */
        final INDArray rot_W_tl_new;
        if (params.isOrthogonalizeAndSortBiasCovariatesEnabled()) {
            orthogonalizeAndSortBiasCovariates(W_tl_new.transpose().mmul(W_tl_new));
            cacheWorkers("after orthogonalization of bias covariates");
            rot_W_tl_new = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl, 0);
        } else {
            rot_W_tl_new = W_tl_new;
        }

        final double errorNormInfinity = CoverageModelEMWorkspaceMathUtils.getINDArrayNormInfinity(
                rot_W_tl_new.sub(W_tl_old));

        /* send the signal to workers for consistency */
        final SubroutineSignal newSig = SubroutineSignal.builder()
                .put("error_norm", errorNormInfinity)
                .put("iterations", iters).build();
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

        switch (params.getWSolverType()) {
            case W_SOLVER_LOCAL:
                /* fetch the required INDArrays */
                Q_ll = mapWorkersAndReduce(cb -> cb.getINDArrayFromCache(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.sum_Q_ll),
                        INDArray::add).div(numTargets);
                Q_tll = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Q_tll, 0);
                Z_ll = sampleBiasLatentPosteriorSecondMoments.sum(0);

                /* instantiate the local implementation of linear operators */
                linop = new CoverageModelWLinearOperatorLocal(Q_tll, Z_ll, regularizerFourierLinearOperator);
                precond = new CoverageModelWPreconditionerLocal(Q_ll, Z_ll, regularizerFourierLinearOperator, numTargets);

                return ImmutablePair.of(linop, precond);

            case W_SOLVER_SPARK:
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
     * Orthogonalize bias covariates and sort the diagonal covariance entries in descending order by performing
     * a rotation in the latent space. This transformation affects:
     *
     *   - worker node copies of E[z_s] and E[z_s z_s^T], and their blocks of W and F[W]
     *   - driver node copies of E[z_s] and E[z_s z_s^T]
     *
     * @param WTW [W]^T [W]
     */
    private void orthogonalizeAndSortBiasCovariates(@Nonnull final INDArray WTW) {
        final ImmutablePair<INDArray, INDArray> ortho =
                CoverageModelEMWorkspaceMathUtils.getOrthogonalizerAndSorterTransformation(WTW, true, logger);
        final INDArray eigs = ortho.left;
        final INDArray U = ortho.right;

        /* update workers */
        pushToWorkers(U, (rot, cb) -> cb.cloneWithRotatedLatentSpace(rot));

        /* update driver node */
        IntStream.range(0, numSamples).parallel().forEach(si -> {
            sampleBiasLatentPosteriorFirstMoments.get(NDArrayIndex.point(si), NDArrayIndex.all())
                    .assign(U.mmul(sampleBiasLatentPosteriorFirstMoments.get(NDArrayIndex.point(si), NDArrayIndex.all())
                            .transpose()).transpose());
            sampleBiasLatentPosteriorSecondMoments.get(NDArrayIndex.point(si), NDArrayIndex.all(), NDArrayIndex.all())
                    .assign(U.mmul(sampleBiasLatentPosteriorSecondMoments.get(NDArrayIndex.point(si), NDArrayIndex.all(), NDArrayIndex.all()))
                            .mmul(U.transpose()));
        });
        biasCovariatesNorm2.assign(eigs);
    }

    /**
     * Same as {@link #orthogonalizeAndSortBiasCovariates(INDArray)} but fetches [W]^T [W] from the workers
     */
    private void orthogonalizeAndSortBiasCovariates() {
        final INDArray WTW = mapWorkersAndReduce(CoverageModelEMComputeBlock::getBiasCovariatesInnerProduct,
                INDArray::add);
        orthogonalizeAndSortBiasCovariates(WTW);
    }

    /**
     * Fetch the log likelihood from compute block(s)
     *
     * @return log likelihood normalized per sample per target
     */
    @EvaluatesRDD
    public double getLogLikelihood() {
        return Arrays.stream(getLogLikelihoodPerSample()).reduce((a, b) -> a + b).orElse(Double.NaN) / numSamples;
    }

    /**
     * Fetch the log likelihood from compute block(s)
     *
     * @return log likelihood normalized per sample per target
     */
    @EvaluatesRDD @CachesRDD
    public double[] getLogLikelihoodPerSample() {
        updateLogLikelihoodCaches();

        final INDArray biasPriorContrib_s = sampleBiasLatentPosteriorFirstMoments
                .mul(sampleBiasLatentPosteriorFirstMoments).muli(-0.5).sum(1);

        final CoverageModelEMComputeBlock.CoverageModelICGCacheNode key;
        if (params.fourierRegularizationEnabled()) {
            key = CoverageModelEMComputeBlock.CoverageModelICGCacheNode.loglike_reg;
        } else {
            key = CoverageModelEMComputeBlock.CoverageModelICGCacheNode.loglike_unreg;
        }
        final INDArray restContrib_s = mapWorkersAndReduce(cb -> cb.getINDArrayFromCache(key), INDArray::add);
        final INDArray sum_M_s = mapWorkersAndReduce(cb -> cb.getINDArrayFromCache(
                CoverageModelEMComputeBlock.CoverageModelICGCacheNode.sum_M_s), INDArray::add);
        final INDArray logLikelihood_s = restContrib_s.add(biasPriorContrib_s).divi(sum_M_s);
        return logLikelihood_s.data().asDouble();
    }

    /**
     * Updates log likelihood caches on compute blocks
     */
    @UpdatesRDD @CachesRDD
    public void updateLogLikelihoodCaches() {
        if (params.fourierRegularizationEnabled()) {
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
                    .map(tb -> new Tuple2<>(tb, new CoverageModelEMComputeBlock(tb, numSamples, numLatents)))
                    .collect(Collectors.toList()), numTargetBlocks)
                    .partitionBy(new HashPartitioner(numTargetBlocks))
                    .cache();
        } else {
            logger.info("Initializing a local compute block");
            localComputeBlock = new CoverageModelEMComputeBlock(targetBlocks.get(0), numSamples, numLatents);
        }
        prevCheckpointedComputeRDD = null;
        cacheCallCounter = 0;
    }

    /**
     * A generic function for handling a blockified list of objects to their corresponding compute nodes
     *
     * If Spark is enabled:
     *
     *      Joins an instance of {@code List<Tuple2<LinearSpaceBlock, V>>} with {@link #computeRDD}, calls the provided
     *      map {@code mapper} on the RDD, and the reference to the old RDD will be replaced with the new RDD.
     *
     * If Spark is disabled:
     *
     *      Only a single target-space block is assumed, such that {@code secondary} is a singleton. The map function
     *      {@code mapper} will be called on the value contained in {@code seconday} and {@link #localComputeBlock}, and
     *      the old instace of {@link CoverageModelEMComputeBlock} is replaced with the new instance returned
     *      by {@code mapper.}
     *
     * @param data the list to joined and mapped together with the compute block(s)
     * @param mapper a mapper binary function that takes a compute block together with an object of type {@code V} and
     *               returns a new compute block
     * @param <V> the type of the object to the broadcasted
     */
    @UpdatesRDD
    private <V> void joinWithWorkersAndMap(@Nonnull final List<Tuple2<LinearSpaceBlock, V>> data,
                                           @Nonnull final Function<Tuple2<CoverageModelEMComputeBlock, V>, CoverageModelEMComputeBlock> mapper) {
        if (sparkContextIsAvailable) {
            final JavaPairRDD<LinearSpaceBlock, V> newRDD =
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
                ex.printStackTrace();
                throw new RuntimeException("Can not apply the map function to the local compute block: " + ex.getMessage());
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
                ex.printStackTrace();
                throw new RuntimeException("Can not apply the map function to the local compute block: " + ex.getMessage());
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
                ex.printStackTrace();
                throw new RuntimeException("Can not apply the map function to the local compute block: " + ex.getMessage());
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
                ex.printStackTrace();
                throw new RuntimeException("Can not apply the map function to the local compute block: " + ex.getMessage());
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
            if (params.isRDDCheckpointingEnabled()) {
                if (cacheCallCounter == params.getRDDCheckpointingInterval()) {
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
    private List<Tuple2<LinearSpaceBlock, INDArray>> chopINDArrayToBlocks(final INDArray arr) {
        if (sparkContextIsAvailable) {
            return CoverageModelSparkUtils.partitionINDArrayToAList(targetBlocks, arr);
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
    private Map<LinearSpaceBlock, INDArray> mapINDArrayToBlocks(final INDArray arr) {
        if (sparkContextIsAvailable) {
            return CoverageModelSparkUtils.partitionINDArrayToAMap(targetBlocks, arr);
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
    private Map<LinearSpaceBlock, ImmutablePair<INDArray, INDArray>> mapINDArrayPairToBlocks(final INDArray arr1,
                                                                                             final INDArray arr2) {
        if (sparkContextIsAvailable) {
            final Map<LinearSpaceBlock, INDArray> map1 =
                    CoverageModelSparkUtils.partitionINDArrayToAMap(targetBlocks, arr1);
            final Map<LinearSpaceBlock, INDArray> map2 =
                    CoverageModelSparkUtils.partitionINDArrayToAMap(targetBlocks, arr2);
            final Map<LinearSpaceBlock, ImmutablePair<INDArray, INDArray>> res = new HashMap<>();
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
     * Returns a {@link Stream<LinearSpaceBlock>} of target-space blocks
     *
     * @return {@link Stream<LinearSpaceBlock>}
     */
    private Stream<LinearSpaceBlock> targetBlockStream() { return targetBlocks.stream(); }

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
        return fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.tot_Psi_st, 1);
    }

    /**
     * Fetches total covariate bias [W.z]_{st} as a sample-target matrix from compute blocks
     *
     * @return an {@link INDArray}
     */
    private INDArray fetchTotalCovariateBiasPerSample() {
        return fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Wz_st, 1);
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
     * Fetches bias covariates from compute blocks
     *
     * @return an {@link INDArray}
     */
    private INDArray fetchBiasCovariates() {
        return fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.W_tl, 0);
    }

    /**
     * Fetches the maximum likelihood estimate of copy ratios as a sample-target matrix from compute blocks
     *
     * @return an {@link INDArray}
     */
    private ImmutablePair<INDArray, INDArray> fetchCopyRatioMaxLikelihoodEstimateData() {

        final INDArray M_Psi_inv_st = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.M_Psi_inv_st, 1);
        final INDArray log_n_st = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.log_n_st, 1);
        final INDArray m_t = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.m_t, 1);
        final INDArray Wz_st = fetchFromWorkers(CoverageModelEMComputeBlock.CoverageModelICGCacheNode.Wz_st, 1);

        /* calculate the required quantities */
        return ImmutablePair.of(log_n_st.sub(Wz_st).subiRowVector(m_t).subiColumnVector(sampleMeanLogReadDepths),
                M_Psi_inv_st);
    }

    /**
     * Fetches the Viterbi copy ratio (or copy number) states as a sample-target matrix
     *
     * @param copyRatioHMMResult an instance of {@link CopyRatioHiddenMarkovModelResults}
     * @param sampleTargetCollections list of target collections for each sample
     * @return an {@link INDArray}
     */
    protected INDArray getViterbiAsNDArray(final List<CopyRatioHiddenMarkovModelResults<CoverageModelCopyRatioEmissionData, S>> copyRatioHMMResult,
                                           final List<TargetCollection<Target>> sampleTargetCollections) {
        final INDArray res = Nd4j.create(numSamples, numTargets);
        for (int ti = 0; ti < numTargets; ti++) {
            final Target target = processedTargetList.get(ti);
            for (int si = 0; si < numSamples; si++) {
                final TargetCollection<Target> sampleTargets = sampleTargetCollections.get(si);
                final List<S> sampleCalls = copyRatioHMMResult.get(si).getViterbiResult();
                final int sampleTargetIndex = sampleTargets.index(target);
                if (sampleTargetIndex >= 0) {
                    res.put(si, ti, sampleCalls.get(sampleTargetIndex).getScalar());
                } else {
                    res.put(si, ti, 0);
                }
            }
        }
        return res;
    }


    /**
     * Saves the model to disk
     *
     * @param outputPath path to write to the model
     */
    public void saveModel(@Nonnull final String outputPath) {
        logger.info("Saving the model to disk...");
        CoverageModelParameters.write(new CoverageModelParameters(processedTargetList,
                fetchMeanLogBias(), fetchTargetUnexplainedVariance(), fetchBiasCovariates()), outputPath);
    }

    /**
     * Saves posteriors to disk
     *
     * @param outputPath path to write posteriors
     * @param verbosityLevel verbosity level
     */
    public void savePosteriors(final String outputPath, final PosteriorVerbosityLevel verbosityLevel) {
        /* create output directory if it doesn't exist */
        createOutputPath(outputPath);

        saveReadDepthPosteriors(outputPath);
        saveLogLikelihoodPosteriors(outputPath);
        saveSampleSpecificUnexplainedVariancePosteriors(outputPath);
        saveBiasLatentPosteriors(outputPath);
        saveTargets(outputPath);

        if (verbosityLevel.equals(PosteriorVerbosityLevel.EXTENDED)) {
            saveCopyRatioPosteriors(outputPath);
            if (params.extendedPosteriorOutputEnabled()) {
                saveCopyRatioMaxLikelihoodEstimates(outputPath);
                saveExtendedPosteriors(outputPath);
            }
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
    protected void saveTargets(final String outputPath) {
        final File targetListFile = new File(outputPath, CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE);
        TargetWriter.writeTargetsToFile(targetListFile, processedTargetList);
    }

    /**
     * Saves copy ratio (or copy number) max likelihood estimates to disk
     *
     * @param outputPath the output path
     */
    protected void saveCopyRatioMaxLikelihoodEstimates(final String outputPath) {
        final ImmutablePair<INDArray, INDArray> copyRatioMLEData = fetchCopyRatioMaxLikelihoodEstimateData();
        final List<String> sampleNames = processedReadCounts.columnNames();
        final List<String> targetNames = processedReadCounts.targets().stream()
                .map(Target::getName).collect(Collectors.toList());

        final File copyRatioMLEFile = new File(outputPath, CoverageModelGlobalConstants.COPY_RATIO_MAX_LIKELIHOOD_ESTIMATES_FILENAME);
        Nd4jIOUtils.writeNDArrayToTextFile(copyRatioMLEData.left, copyRatioMLEFile, sampleNames, targetNames);

        final File copyRatioPrecisionFile = new File(outputPath, CoverageModelGlobalConstants.COPY_RATIO_PRECISION_FILENAME);
        Nd4jIOUtils.writeNDArrayToTextFile(copyRatioMLEData.right, copyRatioPrecisionFile, sampleNames, targetNames);
    }

    /**
     * Saves read depth posteriors to disk
     *
     * @param outputPath the output path
     */
    protected void saveReadDepthPosteriors(final String outputPath) {
        final List<String> sampleNames = processedReadCounts.columnNames();
        final INDArray combinedReadDepthPosteriors = Nd4j.hstack(sampleMeanLogReadDepths, sampleVarLogReadDepths);
        final File sampleReadDepthPosteriorsFile = new File(outputPath, CoverageModelGlobalConstants.SAMPLE_READ_DEPTH_POSTERIORS_FILENAME);
        Nd4jIOUtils.writeNDArrayToTextFile(combinedReadDepthPosteriors, sampleReadDepthPosteriorsFile,
                sampleNames, Arrays.asList("READ_DEPTH_MEAN", "READ_DEPTH_VAR"));
    }

    /**
     * Saves model mog likelihood posteriors to disk
     *
     * @param outputPath the output path
     */
    protected void saveLogLikelihoodPosteriors(final String outputPath) {
        final List<String> sampleNames = processedReadCounts.columnNames();
        final File sampleLogLikelihoodsFile = new File(outputPath, CoverageModelGlobalConstants.SAMPLE_LOG_LIKELIHOODS_FILENAME);
        final INDArray sampleLogLikelihoods = Nd4j.create(getLogLikelihoodPerSample(), new int[] {numSamples, 1});
        Nd4jIOUtils.writeNDArrayToTextFile(sampleLogLikelihoods, sampleLogLikelihoodsFile,
                sampleNames, Collections.singletonList("LOG_LIKELIHOOD"));
    }

    /**
     * Saves sample-specific unexplained variance to disk
     *
     * @param outputPath the output path
     */
    protected void saveSampleSpecificUnexplainedVariancePosteriors(final String outputPath) {
        final List<String> sampleNames = processedReadCounts.columnNames();
        final File sampleUnexplainedVarianceFile = new File(outputPath, CoverageModelGlobalConstants.SAMPLE_UNEXPLAINED_VARIANCE_FILENAME);
        Nd4jIOUtils.writeNDArrayToTextFile(sampleUnexplainedVariance, sampleUnexplainedVarianceFile,
                sampleNames, Collections.singletonList("SAMPLE_UNEXPLAINED_VARIANCE"));
    }

    /**
     * Saves bias latent posteriors E[z_{s\mu}] to disk
     *
     * @param outputPath the output path
     */
    protected void saveBiasLatentPosteriors(final String outputPath) {
        final List<String> sampleNames = processedReadCounts.columnNames();
        final File sampleBiasLatentPosteriorsFile = new File(outputPath,
                CoverageModelGlobalConstants.SAMPLE_BIAS_LATENT_POSTERIORS_FILENAME);
        Nd4jIOUtils.writeNDArrayToTextFile(sampleBiasLatentPosteriorFirstMoments, sampleBiasLatentPosteriorsFile,
                sampleNames, IntStream.range(0, numLatents)
                        .mapToObj(li -> String.format("PC_%d", li)).collect(Collectors.toList()));
    }

    /**
     * Saves copy ratio posteriors (segments, VCF, etc.) to disk
     *
     * @param outputPath the output path
     */
    protected abstract void saveCopyRatioPosteriors(final String outputPath);

    /**
     * Saves extended posteriors to disk
     *
     * @param outputPath the output path
     */
    protected void saveExtendedPosteriors(final String outputPath) {
        final List<String> sampleNames = processedReadCounts.columnNames();
        final List<String> targetNames = processedReadCounts.targets().stream()
                .map(Target::getName).collect(Collectors.toList());

        /* save total unexplained variance as a matrix */
        final File totalExplainedVarianceFile = new File(outputPath, CoverageModelGlobalConstants.TOTAL_UNEXPLAINED_VARIANCE_FILENAME);
        Nd4jIOUtils.writeNDArrayToTextFile(fetchTotalUnexplainedVariance(), totalExplainedVarianceFile,
                sampleNames, targetNames);

        /* save total covariate bias per sample as a matrix */
        final File totalCovariateBiasFile = new File(outputPath, CoverageModelGlobalConstants.TOTAL_COVARIATE_BIAS_FILENAME);
        Nd4jIOUtils.writeNDArrayToTextFile(fetchTotalCovariateBiasPerSample(), totalCovariateBiasFile,
                sampleNames, targetNames);
    }
}
