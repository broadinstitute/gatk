package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Integration test for {@link PostprocessGermlineCNVCalls}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class PostprocessGermlineCNVCallsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/gcnv-postprocess");

    private static final List<String> CALL_SHARDS = Arrays.asList(
            new File(TEST_SUB_DIR, "shard_0-calls").getAbsolutePath(),
            new File(TEST_SUB_DIR, "shard_1-calls").getAbsolutePath(),
            new File(TEST_SUB_DIR, "shard_2-calls").getAbsolutePath());

    private static final List<String> MODEL_SHARDS = Arrays.asList(
            new File(TEST_SUB_DIR, "shard_0-model").getAbsolutePath(),
            new File(TEST_SUB_DIR, "shard_1-model").getAbsolutePath(),
            new File(TEST_SUB_DIR, "shard_2-model").getAbsolutePath());

    private static final String PLOIDY_CALLS = new File(TEST_SUB_DIR, "ploidy-calls").getAbsolutePath();

    private static final int AUTOSOMAL_REF_COPY_NUMBER = 2;

    private final List<String> ALLOSOMAL_CONTIGS = Arrays.asList("X", "Y");

    private static final List<File> INTERVALS_VCF_CORRECT_OUTPUTS = Arrays.asList(
            new File(TEST_SUB_DIR, "intervals_output_SAMPLE_000.vcf"),
            new File(TEST_SUB_DIR, "intervals_output_SAMPLE_001.vcf"),
            new File(TEST_SUB_DIR, "intervals_output_SAMPLE_002.vcf"));

    private static final List<File> SEGMENTS_VCF_CORRECT_OUTPUTS = Arrays.asList(
            new File(TEST_SUB_DIR, "segments_output_SAMPLE_000.vcf"),
            new File(TEST_SUB_DIR, "segments_output_SAMPLE_001.vcf"),
            new File(TEST_SUB_DIR, "segments_output_SAMPLE_002.vcf"));

    private static final int NUM_TEST_SAMPLES = 3;
    private static final int TEST_CALLS_MAX_COPY_NUMBER = 5;

    /**
     * Runs {@link PostprocessGermlineCNVCalls} for a single sample. If {@code segmentsOutputVCF} is null,
     * the tool will only generate intervals VCF output (which is the expected behavior).
     */
    public void runToolForSingleSample(final List<String> callShards,
                                       final List<String> modelShards,
                                       final int sampleIndex,
                                       final File intervalsOutputVCF,
                                       final File segmentsOutputVCF,
                                       final List<String> allosomalContigs,
                                       final int refAutosomalCopyNumber,
                                       final boolean dryRun) {
        final ArgumentsBuilder argumentsBuilder = new ArgumentsBuilder();
        argumentsBuilder.addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE,
                "false");

        argumentsBuilder.addArgument(PostprocessGermlineCNVCalls.SAMPLE_INDEX_LONG_NAME,
                String.valueOf(sampleIndex));
        callShards.forEach(callDir -> argumentsBuilder.addArgument(
                PostprocessGermlineCNVCalls.CALLS_SHARD_PATH_LONG_NAME, callDir));
        argumentsBuilder.addArgument(PostprocessGermlineCNVCalls.OUTPUT_INTERVALS_VCF_LONG_NAME,
                intervalsOutputVCF.getAbsolutePath());
        allosomalContigs.forEach(contig -> argumentsBuilder.addArgument(
                PostprocessGermlineCNVCalls.ALLOSOMAL_CONTIG_LONG_NAME, contig));
        argumentsBuilder.addArgument(PostprocessGermlineCNVCalls.AUTOSOMAL_REF_COPY_NUMBER_LONG_NAME,
                String.valueOf(refAutosomalCopyNumber));

        /* add the rest of the required arguments if segments VCF is required */
        if (segmentsOutputVCF != null) {
            argumentsBuilder.addArgument(PostprocessGermlineCNVCalls.CONTIG_PLOIDY_CALLS_LONG_NAME,
                    PLOIDY_CALLS);
            modelShards.forEach(modelDir -> argumentsBuilder.addArgument(
                    PostprocessGermlineCNVCalls.MODEL_SHARD_PATH_LONG_NAME, modelDir));
            argumentsBuilder.addArgument(PostprocessGermlineCNVCalls.OUTPUT_SEGMENTS_VCF_LONG_NAME,
                    segmentsOutputVCF.getAbsolutePath());
        }
        argumentsBuilder.addArgument("verbosity", "DEBUG");

        if (dryRun) {
            argumentsBuilder.addArgument(PostprocessGermlineCNVCalls.DRY_RUN_LONG_NAME, "true");
        }

        runCommandLine(argumentsBuilder.getArgsList());
    }

    @Test(dataProvider = "differentValidInput", groups = {"python"})
    public void testDifferentValidInputIntervalsComplete(final List<String> callShards,
                                                         final List<String> modelShards) throws IOException {
        for (int sampleIndex = 0; sampleIndex < NUM_TEST_SAMPLES; sampleIndex++) {
            final File actualIntervalsOutputVCF = createTempFile("intervals-output-vcf-" + sampleIndex, ".vcf");
            final File actualSegmentsOutputVCF = createTempFile("segments-output-vcf-" + sampleIndex, ".vcf");
            final File expectedIntervalsOutputVCF = INTERVALS_VCF_CORRECT_OUTPUTS.get(sampleIndex);
            final File expectedSegmentsOutputVCF = SEGMENTS_VCF_CORRECT_OUTPUTS.get(sampleIndex);
            runToolForSingleSample(callShards, modelShards, sampleIndex,
                    actualIntervalsOutputVCF, actualSegmentsOutputVCF,
                    ALLOSOMAL_CONTIGS, AUTOSOMAL_REF_COPY_NUMBER, false);
            IntegrationTestSpec.assertEqualTextFiles(actualIntervalsOutputVCF, expectedIntervalsOutputVCF);
            IntegrationTestSpec.assertEqualTextFiles(actualSegmentsOutputVCF, expectedSegmentsOutputVCF);
        }
    }

    @Test(dataProvider = "differentValidInput")
    public void testDifferentValidInputIntervalsVCFOnly(final List<String> callShards,
                                                        final List<String> modelShards) throws IOException {
        for (int sampleIndex = 0; sampleIndex < NUM_TEST_SAMPLES; sampleIndex++) {
            final File actualIntervalsOutputVCF = createTempFile("intervals-output-vcf-" + sampleIndex, ".vcf");
            final File expectedIntervalsOutputVCF = INTERVALS_VCF_CORRECT_OUTPUTS.get(sampleIndex);
            runToolForSingleSample(callShards, modelShards, sampleIndex,
                    actualIntervalsOutputVCF, null,
                    ALLOSOMAL_CONTIGS, AUTOSOMAL_REF_COPY_NUMBER, false);
            IntegrationTestSpec.assertEqualTextFiles(actualIntervalsOutputVCF, expectedIntervalsOutputVCF);
        }
    }

    @Test(dataProvider = "differentInvalidInput")
    public void testDifferentInvalidInput(final List<String> callShards,
                                          final List<String> modelShards) throws IOException {
        /* dry run -- just input data validation */
        for (int sampleIndex = 0; sampleIndex < NUM_TEST_SAMPLES; sampleIndex++) {
            try {
                runToolForSingleSample(callShards, modelShards, sampleIndex,
                        createTempFile("intervals-output-vcf", ".vcf"),
                        createTempFile("segments-output-vcf", ".vcf"),
                        ALLOSOMAL_CONTIGS, AUTOSOMAL_REF_COPY_NUMBER, true);
            } catch (final IllegalArgumentException ex) {
                continue;
            }
            throw new AssertionError("The test should have failed on invalid input.");
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadAutosomalContigs() {
        runToolForSingleSample(CALL_SHARDS, MODEL_SHARDS, 0,
                createTempFile("intervals-output-vcf", ".vcf"),
                createTempFile("segments-output-vcf", ".vcf"),
                Collections.singletonList("Z"), /* unknown contig */
                AUTOSOMAL_REF_COPY_NUMBER, true);
    }

    @Test(dataProvider = "badRefAutosomalCopyNumber")
    public void testOutOfRangeRefState(final int badRefAutosomalCopyNumber) {
        try {
            runToolForSingleSample(CALL_SHARDS, MODEL_SHARDS, 0,
                    createTempFile("intervals-output-vcf", ".vcf"),
                    createTempFile("segments-output-vcf", ".vcf"),
                    ALLOSOMAL_CONTIGS,
                    badRefAutosomalCopyNumber, /* out of range */
                    true);
        } catch (final CommandLineException.OutOfRangeArgumentValue | IllegalArgumentException ex) {
            return;
        }
        throw new AssertionError("The test should have failed on invalid input.");
    }

    @DataProvider(name = "badRefAutosomalCopyNumber")
    public Object[][] getBadRefAutosomalCopyNumber() {
        return new Object[][] {{-1}, {TEST_CALLS_MAX_COPY_NUMBER + 1}};
    }

    @DataProvider(name = "differentValidInput")
    public Object[][] getDifferentValidInputTestData() {
        return new Object[][] {
                { /* correct order */
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }, { /* permutation */
                        Arrays.asList(CALL_SHARDS.get(2), CALL_SHARDS.get(1), CALL_SHARDS.get(0)),
                        Arrays.asList(MODEL_SHARDS.get(2), MODEL_SHARDS.get(0), MODEL_SHARDS.get(1))
                }, { /* another permutation */
                        Arrays.asList(CALL_SHARDS.get(1), CALL_SHARDS.get(2), CALL_SHARDS.get(0)),
                        Arrays.asList(MODEL_SHARDS.get(1), MODEL_SHARDS.get(0), MODEL_SHARDS.get(2))
                }
        };
    }

    @DataProvider(name = "differentInvalidInput")
    public Object[][] getDifferentInvalidInput() {
        return new Object[][] {
                { /* fewer model shards than call shards */
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1))
                }, { /* fewer call shards than model shards */
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }, { /* repeated calls */
                        Arrays.asList(CALL_SHARDS.get(2), CALL_SHARDS.get(2), CALL_SHARDS.get(1)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }, { /* repeated models */
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(1), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }
        };
    }
}