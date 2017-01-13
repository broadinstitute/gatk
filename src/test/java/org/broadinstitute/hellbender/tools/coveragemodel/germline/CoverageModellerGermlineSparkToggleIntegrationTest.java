package org.broadinstitute.hellbender.tools.coveragemodel.germline;

import htsjdk.samtools.util.Log;
import junit.framework.Assert;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.RealMatrix;
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
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link CoverageModellerGermlineSparkToggle}
 *
 * TODO github/gatk-protected issue #803 test case-sample calling on rearranged targets
 * TODO github/gatk-protected issue #803 test concordance on parameter estimation
 * TODO github/gatk-protected issue #803 test Spark results match local results (of course they do, but is good to test anyway...)
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
    private static final int NUM_LATENTS = 5; /* simulated data uses 3 */
    private static final int MAX_COPY_NUMBER = 3; /* reflects the simulated data */

    private static final int MIN_LEARNING_READ_COUNT = 10;
    private static final int MAX_LEARNING_EM_ITERATIONS = 10;
    private static final int MAX_CALLING_EM_ITERATIONS = 10;

    private static final double MIN_PASS_REF_CONCORDANCE = 0.95;
    private static final double MIN_PASS_ALT_CONCORDANCE = 0.30;
    private static final double MIN_PASS_HOM_DEL_CONCORDANCE = 0.95;
    private static final double MIN_PASS_ABS_CONCORDANCE = 0.30;

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
                "--" + CoverageModelEMParams.NUM_LATENTS_LONG_NAME,
                    String.valueOf(NUM_LATENTS),
                "--" + CoverageModelEMParams.MAPPING_ERROR_RATE_LONG_NAME,
                    String.valueOf(MAPPING_ERROR_RATE),
                "--" + CoverageModelEMParams.MIN_LEARNING_READ_COUNT_LONG_NAME,
                    String.valueOf(MIN_LEARNING_READ_COUNT),
                "--" + CoverageModelEMParams.RUN_CHECKPOINTING_PATH_LONG_NAME,
                    "false",
                "--" + CoverageModelEMParams.NUMBER_OF_TARGET_SPACE_PARTITIONS_LONG_NAME,
                    String.valueOf(SPARK_NUMBER_OF_PARTITIONS),
                "--" + CoverageModellerGermlineSparkToggle.COPY_NUMBER_TRANSITION_PRIOR_TABLE_LONG_NAME,
                    TEST_HMM_PRIORS_TABLE_FILE.getAbsolutePath(),
                "--" + CoverageModellerGermlineSparkToggle.CONTIG_PLOIDY_ANNOTATIONS_TABLE_LONG_NAME,
                    TEST_CONTIG_PLOIDY_ANNOTATIONS_FILE.getAbsolutePath(),
                "--" + CoverageModelEMParams.RDD_CHECKPOINTING_PATH_LONG_NAME,
                    SPARK_CHECKPOINTING_PATH.getAbsolutePath(),
                "--verbosity", "INFO"
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
                    String.valueOf(MAX_LEARNING_EM_ITERATIONS)
        }, getBaseArgs(extraArgs));
    }

    private String[] getCallingArgs(final String... extraArgs) {
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
                    String.valueOf(MAX_CALLING_EM_ITERATIONS)
        }, getBaseArgs(extraArgs));
    }

    private void runLearningAndCallingTest(final String... extraArgs) {
        runCommandLine(getLearningArgs(extraArgs));
        final List<Target> modelledTargets = TargetTableReader.readTargetFile(new File(LEARNING_MODEL_OUTPUT_PATH,
                CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE));

        assertInferredCopyNumberConcordance(LEARNING_POSTERIORS_OUTPUT_PATH, TEST_LEARNING_COMBINED_COPY_NUMBER_FILE,
                modelledTargets, LEARNING_SEX_GENOTYPES_DATA);
        logger.info("Copy number concordance test passed for simultaneous learning and calling");
    }

    private void runCaseSampleCallingTestOnExactModelParams(final String... extraArgs) {
        runCommandLine(getCallingArgs(ArrayUtils.addAll(new String[] {
                "--" + CoverageModellerGermlineSparkToggle.INPUT_MODEL_PATH_LONG_NAME,
                TEST_TRUTH_SIM_MODEL.getAbsolutePath() }, extraArgs)));

        final List<Target> callingTargets = TargetTableReader.readTargetFile(new File(CALLING_POSTERIORS_OUTPUT_PATH,
                CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE));
        assertInferredCopyNumberConcordance(CALLING_POSTERIORS_OUTPUT_PATH, TEST_CALLING_COMBINED_COPY_NUMBER_FILE,
                callingTargets, CALLING_SEX_GENOTYPES_DATA);
        logger.info("Copy number concordance test passed for case sample calling");
    }

    private void runCaseSampleCallingTestOnLearnedModelParams(final String... extraArgs) {
        runCommandLine(getCallingArgs(ArrayUtils.addAll(new String[] {
                "--" + CoverageModellerGermlineSparkToggle.INPUT_MODEL_PATH_LONG_NAME,
                LEARNING_MODEL_OUTPUT_PATH.getAbsolutePath()}, extraArgs)));

        final List<Target> callingTargets = TargetTableReader.readTargetFile(new File(CALLING_POSTERIORS_OUTPUT_PATH,
                CoverageModelGlobalConstants.TARGET_LIST_OUTPUT_FILE));
        assertInferredCopyNumberConcordance(CALLING_POSTERIORS_OUTPUT_PATH, TEST_CALLING_COMBINED_COPY_NUMBER_FILE,
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
    private void assertInferredCopyNumberConcordance(@Nonnull final File posteriorsOutputPath,
                                                     @Nonnull final File truthCopyNumberFile,
                                                     @Nonnull final List<Target> targets,
                                                     @Nonnull final SexGenotypeDataCollection sexGenotypeDataCollection) {
        final ReadCountCollection truthCopyNumberCollection = loadTruthCopyNumberTable(truthCopyNumberFile,
                targets);

        final RealMatrix inferredCopyNumberMatrix = Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(
                Nd4jIOUtils.readNDArrayFromTextFile(new File(posteriorsOutputPath,
                        CoverageModelGlobalConstants.COPY_RATIO_VITERBI_FILENAME)).transpose());
        final ReadCountCollection inferredCopyNumberCollection = new ReadCountCollection(targets,
                truthCopyNumberCollection.columnNames(), inferredCopyNumberMatrix);

        final int numSamples = inferredCopyNumberCollection.columnNames().size();
        final List<String> sampleSexGenotypes = truthCopyNumberCollection.columnNames().stream()
                .map(sampleName -> sexGenotypeDataCollection.getSampleSexGenotypeData(sampleName).getSexGenotype())
                .collect(Collectors.toList());

        final List<SampleCopyNumberConcordanceSummary> sampleConcordanceSummaryList = IntStream.range(0, numSamples)
                .mapToObj(si -> calculateSampleCopyNumberConcordance(truthCopyNumberCollection,
                        inferredCopyNumberCollection, si, sampleSexGenotypes.get(si)))
                .collect(Collectors.toList());

        /* calculate average concordances */
        final double averageRefConcordance = IntStream.range(0, numSamples)
                .mapToDouble(si -> sampleConcordanceSummaryList.get(si).refConcordance)
                .sum() / numSamples;
        final double averageAltConcordance = IntStream.range(0, numSamples)
                .mapToDouble(si -> sampleConcordanceSummaryList.get(si).altConcordance)
                .sum() / numSamples;
        final double[] averageAbsConcordance = IntStream.range(0, numSamples)
                .mapToObj(si -> sampleConcordanceSummaryList.get(si).absConcordance)
                .reduce(new double[MAX_COPY_NUMBER + 1], (a, b) -> {
                    final double[] res = new double[MAX_COPY_NUMBER + 1];
                    for (int i = 0; i <= MAX_COPY_NUMBER; i++) {
                        res[i] = a[i] + b[i];
                    }
                    return res;
                });
        for (int i = 0; i <= MAX_COPY_NUMBER; i++) {
            averageAbsConcordance[i] /= numSamples;
        }
        final double averageOverallAbsConcordance =
                Arrays.stream(averageAbsConcordance).sum() / (MAX_COPY_NUMBER + 1);

        logger.info(String.format("The ref copy number concordance is %f and passes the expected threshold %f",
                averageRefConcordance, MIN_PASS_REF_CONCORDANCE));
        logger.info(String.format("The alt copy number concordance is %f and passes the expected threshold %f",
                averageAltConcordance, MIN_PASS_ALT_CONCORDANCE));
        logger.info(String.format("The homozygous deletion calling concordance is %f and passes the expected threshold %f",
                averageAbsConcordance[0], MIN_PASS_HOM_DEL_CONCORDANCE));
        logger.info(String.format("The average absolute copy number calling concordance %f and passes the expected threshold %f",
                averageOverallAbsConcordance, MIN_PASS_ABS_CONCORDANCE));

        Assert.assertTrue(String.format("The ref copy number concordance %f is lower than expected threshold %f",
                averageRefConcordance, MIN_PASS_REF_CONCORDANCE), averageRefConcordance > MIN_PASS_REF_CONCORDANCE);

        Assert.assertTrue(String.format("The alt copy number concordance %f is lower than expected threshold %f",
                averageAltConcordance, MIN_PASS_ALT_CONCORDANCE), averageAltConcordance > MIN_PASS_ALT_CONCORDANCE);

        Assert.assertTrue(String.format("The homozygous deletion calling concordance %f is lower than expected threshold %f",
                averageAbsConcordance[0], MIN_PASS_HOM_DEL_CONCORDANCE), averageAbsConcordance[0] > MIN_PASS_HOM_DEL_CONCORDANCE);

        Assert.assertTrue(String.format("The average absolute copy number calling concordance %f is lower than expected threshold %f",
                averageOverallAbsConcordance, MIN_PASS_ABS_CONCORDANCE), averageOverallAbsConcordance > MIN_PASS_ABS_CONCORDANCE);
    }

    private SampleCopyNumberConcordanceSummary calculateSampleCopyNumberConcordance(final ReadCountCollection truth,
                                                                                                                             final ReadCountCollection inferred,
                                                                                                                             final int sampleIndex,
                                                                                                                             final String sampleSexGenotype) {
        final List<Target> targets = truth.targets();
        int targetIndex = 0;

        final int[] truthCopyNumberCallsPerSample = Arrays.stream(truth.getColumn(sampleIndex))
                .mapToInt(d -> (int)d).toArray();
        final int[] inferredCopyNumberCallsPerSample = Arrays.stream(inferred.getColumn(sampleIndex))
                .mapToInt(d -> (int)d).toArray();

        final int[] correctCopyNumberCallCounts = new int[MAX_COPY_NUMBER + 1];
        final int[] truthCopyNumberCallCounts = new int[MAX_COPY_NUMBER + 1];
        int truthRefCallCounts = 0, truthAltCallCounts = 0;
        int correctRefCallCounts = 0, correctAltCallCounts = 0;
        for (final Target target : targets) {
            final int refCopyNumber = GERMLINE_PLOIDY_ANNOTATIONS.getTargetGermlinePloidyByGenotype(target,
                    sampleSexGenotype);
            final int truthCopyNumber = truthCopyNumberCallsPerSample[targetIndex];
            final int inferredCopyNumber = inferredCopyNumberCallsPerSample[targetIndex];

            /* ref/alt-blind concordance */
            truthCopyNumberCallCounts[truthCopyNumber]++;
            if (inferredCopyNumber == truthCopyNumber) {
                correctCopyNumberCallCounts[truthCopyNumber]++;
            }

            /* ref/alt-resolved concordance */
            if (truthCopyNumber == refCopyNumber) { /* ref truth */
                truthRefCallCounts++;
                if (inferredCopyNumber == truthCopyNumber) {
                    correctRefCallCounts++;
                }
            } else { /* alt truth */
                truthAltCallCounts++;
                if (inferredCopyNumber == truthCopyNumber) {
                    correctAltCallCounts++;
                }
            }
            targetIndex++;
        }

        final double[] absConcordance = new double[MAX_COPY_NUMBER + 1];
        for (int cn = 0; cn <= MAX_COPY_NUMBER; cn++) {
            if (truthCopyNumberCallCounts[cn] > 0) {
                absConcordance[cn] = ((double) correctCopyNumberCallCounts[cn]) / truthCopyNumberCallCounts[cn];
            } else {
                absConcordance[cn] = 1.0; /* truth has to be 0 too */
            }
        }
        final double refConcordance = truthRefCallCounts > 0 ?
                ((double)correctRefCallCounts) / truthRefCallCounts : 1.0;
        final double altConcordance = truthAltCallCounts > 0 ?
                ((double)correctAltCallCounts) / truthAltCallCounts : 1.0;

        return new SampleCopyNumberConcordanceSummary(refConcordance, altConcordance, absConcordance);
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

    private final class SampleCopyNumberConcordanceSummary {
        /**
         * Fraction of correct ref calls
         */
        final double refConcordance;

        /**
         * Fraction of correct alt calls
         */
        final double altConcordance;

        /**
         * Fraction of correct calls by copy number
         */
        final double[] absConcordance;

        public SampleCopyNumberConcordanceSummary(final double refConcordance,
                                                  final double altConcordance,
                                                  final double[] absConcordance) {
            this.refConcordance = refConcordance;
            this.altConcordance = altConcordance;
            this.absConcordance = absConcordance;
        }

        @Override
        public String toString() {
            return "SampleCopyNumberConcordanceSummary{" +
                    "refConcordance=" + refConcordance +
                    ", altConcordance=" + altConcordance +
                    ", absConcordance=" + Arrays.toString(absConcordance) +
                    '}';
        }
    }
}
