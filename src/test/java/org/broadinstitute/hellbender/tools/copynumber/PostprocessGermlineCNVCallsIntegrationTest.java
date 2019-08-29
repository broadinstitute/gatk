package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Integration test for {@link PostprocessGermlineCNVCalls}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public final class PostprocessGermlineCNVCallsIntegrationTest extends CommandLineProgramTest {
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

    private static final List<File> DENOISED_COPY_RATIOS_OUTPUTS = Arrays.asList(
            new File(TEST_SUB_DIR, "denoised_copy_ratios_SAMPLE_000.tsv"),
            new File(TEST_SUB_DIR, "denoised_copy_ratios_SAMPLE_001.tsv"),
            new File(TEST_SUB_DIR, "denoised_copy_ratios_SAMPLE_002.tsv"));

    private static final int NUM_TEST_SAMPLES = 3;

    /**
     * Runs {@link PostprocessGermlineCNVCalls} for a single sample. If {@code segmentsOutputVCF} is null,
     * the tool will only generate intervals VCF output (which is the expected behavior).
     */
    private void runToolForSingleSample(final List<String> callShards,
                                        final List<String> modelShards,
                                        final int sampleIndex,
                                        final File intervalsOutputVCF,
                                        final File segmentsOutputVCF,
                                        final File denoisedCopyRatiosOutput,
                                        final List<String> allosomalContigs,
                                        final int refAutosomalCopyNumber) {
        final ArgumentsBuilder argumentsBuilder = new ArgumentsBuilder();
        argumentsBuilder.add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE,
                "false");

        argumentsBuilder.add(PostprocessGermlineCNVCalls.SAMPLE_INDEX_LONG_NAME,
                String.valueOf(sampleIndex));
        callShards.forEach(callDir -> argumentsBuilder.add(
                PostprocessGermlineCNVCalls.CALLS_SHARD_PATH_LONG_NAME, callDir));
        argumentsBuilder.add(PostprocessGermlineCNVCalls.OUTPUT_INTERVALS_VCF_LONG_NAME,
                intervalsOutputVCF.getAbsolutePath());
        allosomalContigs.forEach(contig -> argumentsBuilder.add(
                PostprocessGermlineCNVCalls.ALLOSOMAL_CONTIG_LONG_NAME, contig));
        argumentsBuilder.add(PostprocessGermlineCNVCalls.AUTOSOMAL_REF_COPY_NUMBER_LONG_NAME,
                String.valueOf(refAutosomalCopyNumber));

        /* add required arguments for segments VCF creation */
        argumentsBuilder.add(PostprocessGermlineCNVCalls.CONTIG_PLOIDY_CALLS_LONG_NAME,
                PLOIDY_CALLS);
        modelShards.forEach(modelDir -> argumentsBuilder.add(
                PostprocessGermlineCNVCalls.MODEL_SHARD_PATH_LONG_NAME, modelDir));
        argumentsBuilder.add(PostprocessGermlineCNVCalls.OUTPUT_SEGMENTS_VCF_LONG_NAME,
                segmentsOutputVCF.getAbsolutePath());

        /* add denoised copy ratio output file path */
        argumentsBuilder.add(PostprocessGermlineCNVCalls.OUTPUT_DENOISED_COPY_RATIOS_LONG_NAME,
                denoisedCopyRatiosOutput.getAbsolutePath());

        runCommandLine(argumentsBuilder.getArgsList());
    }

    @Test(dataProvider = "differentValidInput", groups = {"python"})
    public void testDifferentValidInput(final int sampleIndex,
                                        final List<String> callShards,
                                        final List<String> modelShards) throws IOException {
        final File actualIntervalsOutputVCF = createTempFile("intervals-output-vcf-" + sampleIndex, ".vcf");
        final File actualSegmentsOutputVCF = createTempFile("segments-output-vcf-" + sampleIndex, ".vcf");
        final File intervalsIndex = new File(actualIntervalsOutputVCF.getAbsolutePath() + FileExtensions.TRIBBLE_INDEX);
        Assert.assertFalse(intervalsIndex.exists());
        final File segmentsIndex = new File(actualIntervalsOutputVCF.getAbsolutePath() + FileExtensions.TRIBBLE_INDEX);
        Assert.assertFalse(segmentsIndex.exists());
        final File actualDenoisedCopyRatiosOutput = createTempFile("denoised-copy-ratios-output-" + sampleIndex, ".tsv");
        final File expectedIntervalsOutputVCF = INTERVALS_VCF_CORRECT_OUTPUTS.get(sampleIndex);
        final File expectedSegmentsOutputVCF = SEGMENTS_VCF_CORRECT_OUTPUTS.get(sampleIndex);
        final File expectedDenoisedCopyRatiosOutput = DENOISED_COPY_RATIOS_OUTPUTS.get(sampleIndex);
        runToolForSingleSample(callShards, modelShards, sampleIndex,
                actualIntervalsOutputVCF, actualSegmentsOutputVCF, actualDenoisedCopyRatiosOutput,
                ALLOSOMAL_CONTIGS, AUTOSOMAL_REF_COPY_NUMBER);

        Assert.assertTrue(intervalsIndex.exists());
        Assert.assertTrue(segmentsIndex.exists());
        IntegrationTestSpec.assertEqualTextFiles(actualIntervalsOutputVCF, expectedIntervalsOutputVCF, "##");
        final VCFHeader intervalsHeader = VariantContextTestUtils.getVCFHeader(actualIntervalsOutputVCF.getAbsolutePath());
        Assert.assertTrue(intervalsHeader.getContigLines().size() > 0);
        IntegrationTestSpec.assertEqualTextFiles(actualSegmentsOutputVCF, expectedSegmentsOutputVCF, "##");
        final VCFHeader segmentsHeader = VariantContextTestUtils.getVCFHeader(actualIntervalsOutputVCF.getAbsolutePath());
        Assert.assertTrue(segmentsHeader.getContigLines().size() > 0);
        IntegrationTestSpec.assertEqualTextFiles(actualDenoisedCopyRatiosOutput, expectedDenoisedCopyRatiosOutput, "##");
    }

    @Test(dataProvider = "differentInvalidInput", expectedExceptions = IllegalArgumentException.class, groups = {"python"})
    public void testDifferentInvalidInput(final int sampleIndex,
                                          final List<String> callShards,
                                          final List<String> modelShards) {
        runToolForSingleSample(callShards, modelShards, sampleIndex,
                createTempFile("intervals-output-vcf", ".vcf"),
                createTempFile("segments-output-vcf", ".vcf"),
                createTempFile("denoised-copy-ratios-output", ".tsv"),
                ALLOSOMAL_CONTIGS, AUTOSOMAL_REF_COPY_NUMBER);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, groups = {"python"})
    public void testBadAutosomalContigs() {
        runToolForSingleSample(CALL_SHARDS, MODEL_SHARDS, 0,
                createTempFile("intervals-output-vcf", ".vcf"),
                createTempFile("segments-output-vcf", ".vcf"),
                createTempFile("denoised-copy-ratios-output", ".tsv"),
                Collections.singletonList("Z"), /* unknown contig */
                AUTOSOMAL_REF_COPY_NUMBER);
    }

    @DataProvider(name = "differentValidInput")
    public Object[][] getDifferentValidInputTestData() {
        final List<Object[]> testCases = new ArrayList<>();
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* correct order */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* permutation */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(2), CALL_SHARDS.get(1), CALL_SHARDS.get(0)),
                        Arrays.asList(MODEL_SHARDS.get(2), MODEL_SHARDS.get(0), MODEL_SHARDS.get(1))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* permutation */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(1), CALL_SHARDS.get(2), CALL_SHARDS.get(0)),
                        Arrays.asList(MODEL_SHARDS.get(1), MODEL_SHARDS.get(0), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        return testCases.toArray(new Object[testCases.size()][]);
    }

    @DataProvider(name = "differentInvalidInput")
    public Object[][] getDifferentInvalidInput() {
        final List<Object[]> testCases = new ArrayList<>();
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* fewer model shards than call shards */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[] {
                        /* fewer call shards than model shards */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[]{
                        /* repeated calls */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(2), CALL_SHARDS.get(2), CALL_SHARDS.get(1)),
                        Arrays.asList(MODEL_SHARDS.get(0), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        IntStream.range(0, NUM_TEST_SAMPLES)
                .mapToObj(sampleIndex ->  new Object[]{
                        /* repeated models */
                        sampleIndex,
                        Arrays.asList(CALL_SHARDS.get(0), CALL_SHARDS.get(1), CALL_SHARDS.get(2)),
                        Arrays.asList(MODEL_SHARDS.get(1), MODEL_SHARDS.get(1), MODEL_SHARDS.get(2))
                }).forEach(testCases::add);
        return testCases.toArray(new Object[testCases.size()][]);
    }
}
