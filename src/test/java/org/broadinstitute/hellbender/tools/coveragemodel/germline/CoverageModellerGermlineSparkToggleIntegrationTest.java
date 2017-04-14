package org.broadinstitute.hellbender.tools.coveragemodel.germline;

import htsjdk.samtools.util.Log;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.AbstractUnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelEMParams;
import org.broadinstitute.hellbender.tools.coveragemodel.CoverageModelGlobalConstants;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jApacheAdapterUtils;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.ContigGermlinePloidyAnnotationTableReader;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.GermlinePloidyAnnotatedTargetCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeDataCollection;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link CoverageModellerGermlineSparkToggle}
 *
 * TODO github/gatk-protected issue #803 test case-sample calling on rearranged targets
 * TODO github/gatk-protected issue #803 test concordance on parameter estimation
 * TODO github/gatk-protected issue #803 test Spark results match local results
 * TODO github/gatk-protected issue #803 report statistics on max likelihood copy ratio as well
 * TODO github/gatk-protected issue #803 report statistics on local copy ratio posteriors as well
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CoverageModellerGermlineSparkToggleIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/coveragemodel";
    private static final File TEST_CONTIG_PLOIDY_ANNOTATIONS_FILE = new File(TEST_SUB_DIR,
            "sim_contig_anots.tsv");
    private static final File TEST_HMM_PRIORS_TABLE_FILE = new File(TEST_SUB_DIR,
            "sim_HMM_priors_table.tsv");

    private static final File TEST_LEARNING_COMBINED_RAW_READ_COUNTS_FILE = new File(TEST_SUB_DIR,
            "learning_combined_read_counts.tsv");
    private static final File TEST_LEARNING_SAMPLE_SEX_GENOTYPES_FILE = new File(TEST_SUB_DIR,
            "learning_sample_sex_genotypes.tsv");
    private static final File TEST_LEARNING_COMBINED_COPY_NUMBER_FILE = new File(TEST_SUB_DIR,
            "learning_combined_copy_number.tsv");

    private static final File TEST_CALLING_COMBINED_RAW_READ_COUNTS_FILE = new File(TEST_SUB_DIR,
            "calling_combined_read_counts.tsv");
    private static final File TEST_CALLING_SAMPLE_SEX_GENOTYPES_FILE = new File(TEST_SUB_DIR,
            "calling_sample_sex_genotypes.tsv");
    private static final File TEST_CALLING_COMBINED_COPY_NUMBER_FILE = new File(TEST_SUB_DIR,
            "calling_combined_copy_number.tsv");

    private static final File TEST_TRUTH_SIM_MODEL = new File(TEST_SUB_DIR, "sim_model");
    private static final File TEST_TARGETS_FILE = new File(TEST_TRUTH_SIM_MODEL, "targets.tsv");

    private static final double MAPPING_ERROR_RATE = 5e-4; /* reflects the simulated data */
    private static final int NUM_TRUTH_LATENTS = 5;

    private static final int NUM_LEARNING_LATENTS = 16;
    private static final boolean ENABLE_LEARNING_GAMMA = false;
    private static final boolean ENABLE_CALLING_GAMMA = false;
    private static final boolean ENABLE_LEARNING_ARD = true;
    private static final CoverageModelEMParams.PsiUpdateMode PSI_UPDATE_MODE =
            CoverageModelEMParams.PsiUpdateMode.PSI_ISOTROPIC;
    private static final CoverageModelEMParams.ModelInitializationStrategy MODEL_INITIALIZATION_STRATEGY =
            CoverageModelEMParams.ModelInitializationStrategy.PCA;
    private static final boolean ENABLE_ADAPTIVE_PSI_UPDATE_SWITCHING = false;

    private static final int MIN_LEARNING_READ_COUNT = 1;
    private static final int MAX_LEARNING_EM_ITERATIONS = 20;
    private static final int MAX_CALLING_EM_ITERATIONS = 10;

    private static final File CHECKPOINTING_PATH = createTempDir("coverage_modeller_germline_checkpointing");
    private static final File LEARNING_OUTPUT_PATH = createTempDir("coverage_modeller_germline_learning_output");
    private static final File CALLING_OUTPUT_PATH = createTempDir("coverage_modeller_germline_calling_output");

    private static final File LEARNING_MODEL_OUTPUT_PATH = new File(LEARNING_OUTPUT_PATH,
            CoverageModellerGermlineSparkToggle.FINAL_MODEL_SUBDIR);
    private static final File LEARNING_POSTERIORS_OUTPUT_PATH = new File(LEARNING_OUTPUT_PATH,
            CoverageModellerGermlineSparkToggle.FINAL_POSTERIORS_SUBDIR);
    private static final File CALLING_POSTERIORS_OUTPUT_PATH = new File(CALLING_OUTPUT_PATH,
            CoverageModellerGermlineSparkToggle.FINAL_POSTERIORS_SUBDIR);

    /* for Spark tests */
    private static final int SPARK_NUMBER_OF_PARTITIONS = 7;
    private static final File SPARK_CHECKPOINTING_PATH = createTempDir("coverage_model_spark_checkpoint");

    private static GermlinePloidyAnnotatedTargetCollection GERMLINE_PLOIDY_ANNOTATIONS;
    private static SexGenotypeDataCollection LEARNING_SEX_GENOTYPES_DATA, CALLING_SEX_GENOTYPES_DATA;

    @BeforeSuite @Override
    public void setTestVerbosity(){
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
    }

    @BeforeSuite
    public void init() throws IOException {
        LEARNING_SEX_GENOTYPES_DATA = new SexGenotypeDataCollection(TEST_LEARNING_SAMPLE_SEX_GENOTYPES_FILE);
        CALLING_SEX_GENOTYPES_DATA = new SexGenotypeDataCollection(TEST_CALLING_SAMPLE_SEX_GENOTYPES_FILE);
        GERMLINE_PLOIDY_ANNOTATIONS = new GermlinePloidyAnnotatedTargetCollection(ContigGermlinePloidyAnnotationTableReader
                .readContigGermlinePloidyAnnotationsFromFile(TEST_CONTIG_PLOIDY_ANNOTATIONS_FILE),
                TargetTableReader.readTargetFile(TEST_TARGETS_FILE));
    }

    private String[] getBaseArgs(final String... extraArgs) {
        return ArrayUtils.addAll(new String[] {
                "--" + CoverageModelEMParams.MAPPING_ERROR_RATE_LONG_NAME,
                    String.valueOf(MAPPING_ERROR_RATE),
                "--" + CoverageModelEMParams.RUN_CHECKPOINTING_PATH_LONG_NAME,
                    CHECKPOINTING_PATH.getAbsolutePath(),
                "--" + CoverageModelEMParams.NUMBER_OF_TARGET_SPACE_PARTITIONS_LONG_NAME,
                    String.valueOf(SPARK_NUMBER_OF_PARTITIONS),
                "--" + CoverageModellerGermlineSparkToggle.COPY_NUMBER_TRANSITION_PRIOR_TABLE_LONG_NAME,
                    TEST_HMM_PRIORS_TABLE_FILE.getAbsolutePath(),
                "--" + CoverageModellerGermlineSparkToggle.CONTIG_PLOIDY_ANNOTATIONS_TABLE_LONG_NAME,
                    TEST_CONTIG_PLOIDY_ANNOTATIONS_FILE.getAbsolutePath(),
                "--" + CoverageModelEMParams.RDD_CHECKPOINTING_PATH_LONG_NAME,
                    SPARK_CHECKPOINTING_PATH.getAbsolutePath(),
                "--" + CoverageModelEMParams.EXTENDED_POSTERIOR_OUTPUT_ENABLED_LONG_NAME,
                    "true",
                "--verbosity",
                    "INFO"
        }, extraArgs);
    }

    private String[] getLearningArgs(final String... extraArgs) {
        return ArrayUtils.addAll(new String[] {
                "--" + CoverageModellerGermlineSparkToggle.JOB_TYPE_LONG_NAME,
                    CoverageModellerGermlineSparkToggle.JobType.LEARN_AND_CALL.name(),
                "--" + CoverageModellerGermlineSparkToggle.INPUT_READ_COUNTS_TABLE_LONG_NAME,
                    TEST_LEARNING_COMBINED_RAW_READ_COUNTS_FILE.getAbsolutePath(),
                "--" + CoverageModellerGermlineSparkToggle.SAMPLE_SEX_GENOTYPE_TABLE_LONG_NAME,
                    TEST_LEARNING_SAMPLE_SEX_GENOTYPES_FILE.getAbsolutePath(),
                "--" + CoverageModellerGermlineSparkToggle.OUTPUT_PATH_LONG_NAME,
                    LEARNING_OUTPUT_PATH.getAbsolutePath(),
                "--" + CoverageModelEMParams.MAX_EM_ITERATIONS_LONG_NAME,
                    String.valueOf(MAX_LEARNING_EM_ITERATIONS),
                "--" + CoverageModelEMParams.MIN_LEARNING_READ_COUNT_LONG_NAME,
                    String.valueOf(MIN_LEARNING_READ_COUNT),
                "--" + CoverageModelEMParams.RUN_CHECKPOINTING_ENABLED_LONG_NAME,
                    "false",
                "--" + CoverageModelEMParams.GAMMA_UPDATE_ENABLED_LONG_NAME,
                    String.valueOf(ENABLE_LEARNING_GAMMA),
                "--" + CoverageModelEMParams.ARD_ENABLED_LONG_NAME,
                    String.valueOf(ENABLE_LEARNING_ARD),
                "--" + CoverageModelEMParams.NUM_LATENTS_LONG_NAME,
                    String.valueOf(NUM_LEARNING_LATENTS),
                "--" + CoverageModelEMParams.ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_LONG_NAME,
                    String.valueOf(ENABLE_ADAPTIVE_PSI_UPDATE_SWITCHING),
                "--" + CoverageModelEMParams.PSI_SOLVER_MODE_LONG_NAME,
                    PSI_UPDATE_MODE.name(),
                "--" + CoverageModelEMParams.MODEL_INITIALIZATION_STRATEGY_LONG_NAME,
                    MODEL_INITIALIZATION_STRATEGY.name()
        }, getBaseArgs(extraArgs));
    }

    private String[] getCallingOnLearnedModelArgs(final String... extraArgs) {
        return ArrayUtils.addAll(new String[] {
                "--" + CoverageModellerGermlineSparkToggle.JOB_TYPE_LONG_NAME,
                    CoverageModellerGermlineSparkToggle.JobType.CALL_ONLY.name(),
                "--" + CoverageModellerGermlineSparkToggle.INPUT_READ_COUNTS_TABLE_LONG_NAME,
                    TEST_CALLING_COMBINED_RAW_READ_COUNTS_FILE.getAbsolutePath(),
                "--" + CoverageModellerGermlineSparkToggle.SAMPLE_SEX_GENOTYPE_TABLE_LONG_NAME,
                    TEST_CALLING_SAMPLE_SEX_GENOTYPES_FILE.getAbsolutePath(),
                "--" + CoverageModellerGermlineSparkToggle.OUTPUT_PATH_LONG_NAME,
                    CALLING_OUTPUT_PATH.getAbsolutePath(),
                "--" + CoverageModelEMParams.MAX_EM_ITERATIONS_LONG_NAME,
                    String.valueOf(MAX_CALLING_EM_ITERATIONS),
                "--" + CoverageModelEMParams.RUN_CHECKPOINTING_ENABLED_LONG_NAME,
                    "false",
                "--" + CoverageModelEMParams.GAMMA_UPDATE_ENABLED_LONG_NAME,
                    String.valueOf(ENABLE_CALLING_GAMMA),
                "--" + CoverageModelEMParams.ARD_ENABLED_LONG_NAME,
                    String.valueOf(ENABLE_LEARNING_ARD),
                "--" + CoverageModelEMParams.NUM_LATENTS_LONG_NAME,
                    String.valueOf(NUM_LEARNING_LATENTS),
                "--" + CoverageModelEMParams.ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_LONG_NAME,
                    String.valueOf(ENABLE_ADAPTIVE_PSI_UPDATE_SWITCHING),
                "--" + CoverageModelEMParams.PSI_SOLVER_MODE_LONG_NAME,
                    PSI_UPDATE_MODE.name()
        }, getBaseArgs(extraArgs));
    }

    private String[] getCallingOnExactModelArgs(final String... extraArgs) {
        return ArrayUtils.addAll(new String[] {
                "--" + CoverageModellerGermlineSparkToggle.JOB_TYPE_LONG_NAME,
                    CoverageModellerGermlineSparkToggle.JobType.CALL_ONLY.name(),
                "--" + CoverageModellerGermlineSparkToggle.INPUT_READ_COUNTS_TABLE_LONG_NAME,
                    TEST_CALLING_COMBINED_RAW_READ_COUNTS_FILE.getAbsolutePath(),
                "--" + CoverageModellerGermlineSparkToggle.SAMPLE_SEX_GENOTYPE_TABLE_LONG_NAME,
                    TEST_CALLING_SAMPLE_SEX_GENOTYPES_FILE.getAbsolutePath(),
                "--" + CoverageModellerGermlineSparkToggle.OUTPUT_PATH_LONG_NAME,
                    CALLING_OUTPUT_PATH.getAbsolutePath(),
                "--" + CoverageModelEMParams.MAX_EM_ITERATIONS_LONG_NAME,
                   String.valueOf(MAX_CALLING_EM_ITERATIONS),
                "--" + CoverageModelEMParams.RUN_CHECKPOINTING_ENABLED_LONG_NAME,
                    "false",
                "--" + CoverageModelEMParams.GAMMA_UPDATE_ENABLED_LONG_NAME,
                    String.valueOf(ENABLE_CALLING_GAMMA),
                "--" + CoverageModelEMParams.ARD_ENABLED_LONG_NAME,
                    "false",
                "--" + CoverageModelEMParams.NUM_LATENTS_LONG_NAME,
                   String.valueOf(NUM_TRUTH_LATENTS),
                "--" + CoverageModelEMParams.ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_LONG_NAME,
                    String.valueOf(ENABLE_ADAPTIVE_PSI_UPDATE_SWITCHING),
                "--" + CoverageModelEMParams.PSI_SOLVER_MODE_LONG_NAME,
                    PSI_UPDATE_MODE.name()
        }, getBaseArgs(extraArgs));
    }

    private void runLearningAndCallingTest(final String... extraArgs) {
        runCommandLine(getLearningArgs(extraArgs));
        final List<Target> modelledTargets = TargetTableReader.readTargetFile(new File(LEARNING_MODEL_OUTPUT_PATH,
                CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE));

        reportCopyNumberSummaryStatistics(LEARNING_POSTERIORS_OUTPUT_PATH, TEST_LEARNING_COMBINED_COPY_NUMBER_FILE,
                modelledTargets, LEARNING_SEX_GENOTYPES_DATA);
        logger.info("Copy number concordance test passed for simultaneous learning and calling");
    }

    private void runCaseSampleCallingTestOnExactModelParams(final String... extraArgs) {
        runCommandLine(getCallingOnExactModelArgs(ArrayUtils.addAll(new String[] {
                "--" + CoverageModellerGermlineSparkToggle.INPUT_MODEL_PATH_LONG_NAME,
                TEST_TRUTH_SIM_MODEL.getAbsolutePath() }, extraArgs)));

        final List<Target> callingTargets = TargetTableReader.readTargetFile(new File(CALLING_POSTERIORS_OUTPUT_PATH,
                CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE));
        reportCopyNumberSummaryStatistics(CALLING_POSTERIORS_OUTPUT_PATH, TEST_CALLING_COMBINED_COPY_NUMBER_FILE,
                callingTargets, CALLING_SEX_GENOTYPES_DATA);
        logger.info("Copy number concordance test passed for case sample calling");
    }

    private void runCaseSampleCallingTestOnLearnedModelParams(final String... extraArgs) {
        runCommandLine(getCallingOnLearnedModelArgs(ArrayUtils.addAll(new String[] {
                "--" + CoverageModellerGermlineSparkToggle.INPUT_MODEL_PATH_LONG_NAME,
                LEARNING_MODEL_OUTPUT_PATH.getAbsolutePath()}, extraArgs)));

        final List<Target> callingTargets = TargetTableReader.readTargetFile(new File(CALLING_POSTERIORS_OUTPUT_PATH,
                CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE));
        reportCopyNumberSummaryStatistics(CALLING_POSTERIORS_OUTPUT_PATH, TEST_CALLING_COMBINED_COPY_NUMBER_FILE,
                callingTargets, CALLING_SEX_GENOTYPES_DATA);
        logger.info("Copy number concordance test passed for case sample calling");
    }

    @Test
    public void runLearningAndCallingTestLocal() {
        runLearningAndCallingTest("--" + SparkToggleCommandLineProgram.DISABLE_SPARK_FULL_NAME, "true");
    }

    @Test(dependsOnMethods = "runLearningAndCallingTestLocal")
    public void runCaseSampleCallingTestOnLearnedModelParamsLocal() {
        runCaseSampleCallingTestOnLearnedModelParams("--" + SparkToggleCommandLineProgram.DISABLE_SPARK_FULL_NAME, "true");
    }

    @Test
    public void runCaseSampleCallingTestOnExactModelParamsLocal() {
        runCaseSampleCallingTestOnExactModelParams("--" + SparkToggleCommandLineProgram.DISABLE_SPARK_FULL_NAME, "true");
    }

    @Test
    public void runLearningAndCallingTestSpark() {
        runLearningAndCallingTest();
    }

    @Test(dependsOnMethods = "runLearningAndCallingTestSpark")
    public void runCaseSampleCallingTestOnLearnedModelParamsSpark() {
        runCaseSampleCallingTestOnLearnedModelParams();
    }

    @Test
    public void runCaseSampleCallingTestOnExactModelParamsSpark() {
        runCaseSampleCallingTestOnExactModelParams();
    }

    /* Shame on me for using {@link ReadCountCollection} to store copy numbers! */
    private void reportCopyNumberSummaryStatistics(@Nonnull final File posteriorsOutputPath,
                                                   @Nonnull final File truthCopyNumberFile,
                                                   @Nonnull final List<Target> targets,
                                                   @Nonnull final SexGenotypeDataCollection sexGenotypeDataCollection) {
        final ReadCountCollection truthCopyNumberCollection = loadTruthCopyNumberTable(truthCopyNumberFile,
                targets);

        final RealMatrix calledCopyNumberMatrix = Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(
                Nd4jIOUtils.readNDArrayMatrixFromTextFile(new File(posteriorsOutputPath,
                        CoverageModelGlobalConstants.COPY_RATIO_VITERBI_FILENAME)));
        final ReadCountCollection calledCopyNumberCollection = new ReadCountCollection(targets,
                truthCopyNumberCollection.columnNames(), calledCopyNumberMatrix);

        final int numSamples = calledCopyNumberCollection.columnNames().size();
        final List<String> sampleSexGenotypes = truthCopyNumberCollection.columnNames().stream()
                .map(sampleName -> sexGenotypeDataCollection.getSampleSexGenotypeData(sampleName).getSexGenotype())
                .collect(Collectors.toList());

        final List<SampleCopyNumberSummaryStatistics> sampleSummaryStatisticsList = IntStream.range(0, numSamples)
                .mapToObj(si -> calculateSampleCopyNumberConcordance(truthCopyNumberCollection,
                        calledCopyNumberCollection, si, sampleSexGenotypes.get(si)))
                .collect(Collectors.toList());

        /* calculation various summary statistics */
        final AbstractUnivariateStatistic calculator = new Mean();
        final ConfusionRates homDelMedianRates = ConfusionMatrix.getConfusionRates(sampleSummaryStatisticsList.stream()
                .map(ss -> ss.homozygousDeletionConfusionMatrix).collect(Collectors.toList()), calculator);
        final ConfusionRates hetDelMedianRates = ConfusionMatrix.getConfusionRates(sampleSummaryStatisticsList.stream()
                .map(ss -> ss.heterozygousDeletionConfusionMatrix).collect(Collectors.toList()), calculator);
        final ConfusionRates dupMedianRates = ConfusionMatrix.getConfusionRates(sampleSummaryStatisticsList.stream()
                .map(ss -> ss.duplicationConfusionMatrix).collect(Collectors.toList()), calculator);
        final double absoluteConcordance = Concordance.getCollectionConcordance(sampleSummaryStatisticsList.stream()
                .map(ss -> ss.absoluteCopyNumberConcordance)
                .collect(Collectors.toList()), calculator);

        /* log */
        logger.info("Homozygous deletion statistics: " + homDelMedianRates.toString());
        logger.info("Heterozygous deletion statistics: " + hetDelMedianRates.toString());
        logger.info("Duplication statistics: " + dupMedianRates.toString());
        logger.info(String.format("Absolute copy number calling concordance: %f", absoluteConcordance));
    }

    private SampleCopyNumberSummaryStatistics calculateSampleCopyNumberConcordance(final ReadCountCollection truth,
                                                                                   final ReadCountCollection called,
                                                                                   final int sampleIndex,
                                                                                   final String sampleSexGenotype) {
        final List<Target> targets = truth.targets();
        final int[] truthCopyNumberCallsPerSample = Arrays.stream(truth.getColumn(sampleIndex))
                .mapToInt(d -> (int)d).toArray();
        final int[] calledCopyNumberCallsPerSample = Arrays.stream(called.getColumn(sampleIndex))
                .mapToInt(d -> (int)d).toArray();
        final int[] refCopyNumberCallsPerSample = targets.stream()
                .mapToInt(target -> GERMLINE_PLOIDY_ANNOTATIONS.getTargetGermlinePloidyByGenotype(target,
                        sampleSexGenotype))
                .toArray();

        return new SampleCopyNumberSummaryStatistics(refCopyNumberCallsPerSample, truthCopyNumberCallsPerSample,
                calledCopyNumberCallsPerSample);
    }

    private ReadCountCollection loadTruthCopyNumberTable(@Nonnull final File truthCopyNumberFile,
                                                         @Nonnull final List<Target> modelledTargets) {
        final ReadCountCollection truthCopyNumberCollection;
        try {
            truthCopyNumberCollection = ReadCountCollectionUtils.parse(truthCopyNumberFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not load truth copy number file");
        }
        return truthCopyNumberCollection.subsetTargets(new LinkedHashSet<>(modelledTargets))
                .arrangeTargets(modelledTargets);
    }

    /**
     * Confusion rates for a binary classification problem
     */
    private static final class ConfusionRates {
        final double TPR, TNR, FPR, FNR;

        ConfusionRates(final double TPR, final double TNR, final double FPR, final double FNR) {
            this.TPR = TPR;
            this.TNR = TNR;
            this.FPR = FPR;
            this.FNR = FNR;
        }

        @Override
        public String toString() {
            return String.format("TPR = %f, FPR = %f, TNR = %f, FNR = %f", TPR, FPR, TNR, FNR);
        }
    }

    /**
     * Confusion matrix of a binary classification problem
     */
    private static final class ConfusionMatrix {
        final int TP, TN, FP, FN;

        /* condition positive = TP + FN */
        final int CP;

        /* condition negative = FP + TN */
        final int CN;

        public static final double DEFAULT_UNDEFINED = -1.0;

        public ConfusionMatrix(final int TP, final int FP, final int TN, final int FN) {
            this.TP = ParamUtils.isPositiveOrZero(TP, "true positive count must be non-negative");
            this.FP = ParamUtils.isPositiveOrZero(FP, "false positive count must be non-negative");
            this.TN = ParamUtils.isPositiveOrZero(TN, "true negative count must be non-negative");
            this.FN = ParamUtils.isPositiveOrZero(FN, "false negative count must be non-negative");
            CP = TP + FN;
            CN = FP + TN;
        }

        public double getTruePositiveRate() {
            return CP > 0 ? ((double) TP) / CP : DEFAULT_UNDEFINED;
        }

        public double getTrueNegativeRate() {
            return CN > 0 ? ((double) TN) / CN : DEFAULT_UNDEFINED;
        }

        public double getFalsePositiveRate() {
            return CN > 0 ? ((double) FP) / CN : DEFAULT_UNDEFINED;
        }

        public double getFalseNegativeRate() {
            return CP > 0 ? ((double) FN) / CP : DEFAULT_UNDEFINED;
        }

        public static double getCollectionTruePositiveRate(@Nonnull final Collection<ConfusionMatrix> col,
                                                           @Nonnull final AbstractUnivariateStatistic calculator) {
            try {
                return calculator.evaluate(col.stream()
                        .mapToDouble(ConfusionMatrix::getTruePositiveRate)
                        .filter(val -> val != ConfusionMatrix.DEFAULT_UNDEFINED)
                        .toArray());
            } catch (final MathIllegalArgumentException ex) {
                return DEFAULT_UNDEFINED;
            }
        }

        public static double getCollectionTrueNegativeRate(@Nonnull final Collection<ConfusionMatrix> col,
                                                           @Nonnull final AbstractUnivariateStatistic calculator) {
            try {
                return calculator.evaluate(col.stream()
                        .mapToDouble(ConfusionMatrix::getTrueNegativeRate)
                        .filter(val -> val != ConfusionMatrix.DEFAULT_UNDEFINED)
                        .toArray());
            } catch (final MathIllegalArgumentException ex) {
                return DEFAULT_UNDEFINED;
            }
        }

        public static double getCollectionFalsePositiveRate(@Nonnull final Collection<ConfusionMatrix> col,
                                                            @Nonnull final AbstractUnivariateStatistic calculator) {
            try {
                return calculator.evaluate(col.stream()
                        .mapToDouble(ConfusionMatrix::getFalsePositiveRate)
                        .filter(val -> val != ConfusionMatrix.DEFAULT_UNDEFINED)
                        .toArray());
            } catch (final MathIllegalArgumentException ex) {
                return DEFAULT_UNDEFINED;
            }
        }

        public static double getCollectionFalseNegativeRate(@Nonnull final Collection<ConfusionMatrix> col,
                                                            @Nonnull final AbstractUnivariateStatistic calculator) {
            try {
                return calculator.evaluate(col.stream()
                        .mapToDouble(ConfusionMatrix::getFalseNegativeRate)
                        .filter(val -> val != ConfusionMatrix.DEFAULT_UNDEFINED)
                        .toArray());
            } catch (final MathIllegalArgumentException ex) {
                return DEFAULT_UNDEFINED;
            }
        }

        public static ConfusionRates getConfusionRates(@Nonnull final Collection<ConfusionMatrix> col,
                                                       @Nonnull final AbstractUnivariateStatistic calculator) {
            return new ConfusionRates(getCollectionTruePositiveRate(col, calculator),
                    getCollectionTrueNegativeRate(col, calculator),
                    getCollectionFalsePositiveRate(col, calculator),
                    getCollectionFalseNegativeRate(col, calculator));
        }

        @Override
        public String toString() {
            return "ConfusionMatrix{" +
                    "TP=" + TP +
                    ", TN=" + TN +
                    ", FP=" + FP +
                    ", FN=" + FN +
                    '}';
        }
    }

    /**
     * Concordance statistics
     */
    private static final class Concordance {
        final int right, wrong;

        public static final double DEFAULT_UNDEFINED = -1.0;

        public Concordance(final int right, final int wrong) {
            this.right = ParamUtils.isPositiveOrZero(right, "Right count must be non-negative");
            this.wrong = ParamUtils.isPositiveOrZero(wrong, "Wrong count must be non-negative");
        }

        public double getConcordance() {
            return (right + wrong > 0) ? (double) right / (right + wrong) : DEFAULT_UNDEFINED;
        }

        public double getDiscordance() {
            return (right + wrong > 0) ? (double) wrong / (right + wrong) : DEFAULT_UNDEFINED;
        }

        public static double getCollectionConcordance(@Nonnull final Collection<Concordance> col,
                                                      @Nonnull final AbstractUnivariateStatistic calculator) {
            try {
                return calculator.evaluate(col.stream()
                        .mapToDouble(Concordance::getConcordance)
                        .filter(val -> val != DEFAULT_UNDEFINED)
                        .toArray());
            } catch (final MathIllegalArgumentException ex) {
                return DEFAULT_UNDEFINED;
            }
        }

        public static double getCollectionDiscordance(@Nonnull final Collection<Concordance> col,
                                                      @Nonnull final AbstractUnivariateStatistic calculator) {
            try {
                return calculator.evaluate(col.stream()
                        .mapToDouble(Concordance::getDiscordance)
                        .filter(val -> val != DEFAULT_UNDEFINED)
                        .toArray());
            } catch (final MathIllegalArgumentException ex) {
                return DEFAULT_UNDEFINED;
            }
        }
    }

    private final class SampleCopyNumberSummaryStatistics {

        final int length;
        final ConfusionMatrix homozygousDeletionConfusionMatrix;
        final ConfusionMatrix heterozygousDeletionConfusionMatrix;
        final ConfusionMatrix duplicationConfusionMatrix;
        final Concordance absoluteCopyNumberConcordance;

        /**
         *
         * @param refCopyNumberArray array of reference copy numbers
         * @param truthCopyNumberArray array of truth copy numbers
         * @param calledCopyNumberArray array of called copy numbers
         */
        SampleCopyNumberSummaryStatistics(final int[] refCopyNumberArray,
                                          final int[] truthCopyNumberArray,
                                          final int[] calledCopyNumberArray) {
            Utils.validateArg(refCopyNumberArray.length == truthCopyNumberArray.length &&
                refCopyNumberArray.length == calledCopyNumberArray.length, "The reference, truth, and called copy number" +
                    " arrays must have the same length");
            length = refCopyNumberArray.length;
            homozygousDeletionConfusionMatrix = calculateHomozygousDeletionStatistics(refCopyNumberArray,
                    truthCopyNumberArray, calledCopyNumberArray);
            heterozygousDeletionConfusionMatrix = calculateHeterozygousDeletionStatistics(refCopyNumberArray,
                    truthCopyNumberArray, calledCopyNumberArray);
            duplicationConfusionMatrix = calculateDuplicationStatistics(refCopyNumberArray,
                    truthCopyNumberArray, calledCopyNumberArray);
            absoluteCopyNumberConcordance = calculateAbsoluteCopyNumberStatistics(truthCopyNumberArray,
                    calledCopyNumberArray);
        }

        private ConfusionMatrix calculateHomozygousDeletionStatistics(final int[] refCopyNumberArray,
                                                                      final int[] truthCopyNumberArray,
                                                                      final int[] calledCopyNumberArray) {
            int TP = 0, FP = 0, TN = 0, FN = 0;
            for (int i = 0; i < length; i++) {
                if (refCopyNumberArray[i] > 0 && truthCopyNumberArray[i] == 0 && calledCopyNumberArray[i] == 0) {
                    TP++;
                }
                if (refCopyNumberArray[i] > 0 && truthCopyNumberArray[i] > 0 && calledCopyNumberArray[i] == 0) {
                    FP++;
                }
                if (truthCopyNumberArray[i] > 0 && calledCopyNumberArray[i] > 0) {
                    TN++;
                }
                if (refCopyNumberArray[i] > 0 && truthCopyNumberArray[i] == 0 && calledCopyNumberArray[i] > 0) {
                    FN++;
                }
            }
            return new ConfusionMatrix(TP, FP, TN, FN);
        }

        private ConfusionMatrix calculateHeterozygousDeletionStatistics(final int[] refCopyNumberArray,
                                                                        final int[] truthCopyNumberArray,
                                                                        final int[] calledCopyNumberArray) {
            int TP = 0, FP = 0, TN = 0, FN = 0;
            for (int i = 0; i < length; i++) {
                if (refCopyNumberArray[i] != 2) {
                    continue;
                }
                if (truthCopyNumberArray[i] == 1 && calledCopyNumberArray[i] == 1) {
                    TP++;
                }
                if (truthCopyNumberArray[i] != 1 && calledCopyNumberArray[i] == 1) {
                    FP++;
                }
                if (truthCopyNumberArray[i] != 1 && calledCopyNumberArray[i] != 1) {
                    TN++;
                }
                if (truthCopyNumberArray[i] == 1 && calledCopyNumberArray[i] != 1) {
                    FN++;
                }
            }
            return new ConfusionMatrix(TP, FP, TN, FN);
        }

        private ConfusionMatrix calculateDuplicationStatistics(final int[] refCopyNumberArray,
                                                               final int[] truthCopyNumberArray,
                                                               final int[] calledCopyNumberArray) {
            int TP = 0, FP = 0, TN = 0, FN = 0;
            for (int i = 0; i < length; i++) {
                if (refCopyNumberArray[i] == 0) {
                    continue;
                }
                if (truthCopyNumberArray[i] > refCopyNumberArray[i] && calledCopyNumberArray[i] > refCopyNumberArray[i]) {
                    TP++;
                }
                if (truthCopyNumberArray[i] <= refCopyNumberArray[i] && calledCopyNumberArray[i] > refCopyNumberArray[i]) {
                    FP++;
                }
                if (truthCopyNumberArray[i] <= refCopyNumberArray[i] && calledCopyNumberArray[i] <= refCopyNumberArray[i]) {
                    TN++;
                }
                if (truthCopyNumberArray[i] > refCopyNumberArray[i] && calledCopyNumberArray[i] <= refCopyNumberArray[i]) {
                    FN++;
                }
            }
            return new ConfusionMatrix(TP, FP, TN, FN);
        }

        private Concordance calculateAbsoluteCopyNumberStatistics(final int[] truthCopyNumberArray,
                                                                  final int[] calledCopyNumberArray) {
            int right = 0, wrong = 0;
            for (int i = 0; i < length; i++) {
                if (truthCopyNumberArray[i] == calledCopyNumberArray[i]) {
                    right++;
                } else {
                    wrong++;
                }
            }
            return new Concordance(right, wrong);
        }
    }

}
