package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.broadinstitute.hellbender.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.PostProcessGermlineCNVCalls;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

/**
 * Integration test for {@link PostProcessGermlineCNVCalls}
 */
public class PostProcessGermlineCNVCallsIntegrationTest extends CommandLineProgramTest {

    // Test directory
    private static final File TEST_DIR = new File(toolsTestDir, "copynumber/gcnv/");

    private static final String SAMPLE_NAME = "SAMPLE_0";

    private static final File VCF_CORRECT_OUTPUT = new File(TEST_DIR, "output.vcf");

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testWrongInputChunkOrder() {
        final List<String> chunkDirectoriesNames = new ArrayList<>(Arrays.asList(
                "/normal_cohort_chunk_2-calls/", "/normal_cohort_chunk_0-calls/", "/normal_cohort_chunk_1-calls/"));
        final List<String> chunkDirectoriesFileList = chunkDirectoriesNames.stream()
                .map(s -> TEST_DIR.getAbsolutePath().concat(s)).collect(Collectors.toList());
        final File outputVCF = createTempFile("test", ".vcf");

        final List<String> arguments = new ArrayList<>();
        chunkDirectoriesFileList.forEach(dir -> {
            arguments.add("-" + PostProcessGermlineCNVCalls.CHUNK_PATH_LONG_NAME); arguments.add(dir);
        });
        arguments.add("-" + PostProcessGermlineCNVCalls.SAMPLE_DIRECTORY_LONG_NAME);
        arguments.add(SAMPLE_NAME);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        arguments.add(outputVCF.getAbsolutePath());

        runCommandLine(arguments);
    }

    @Test
    public void testCorrectGenerateVCFRun() throws IOException {
        final List<String> chunkDirectoriesNames = new ArrayList<>(Arrays.asList(
                "/normal_cohort_chunk_0-calls/", "/normal_cohort_chunk_1-calls/", "/normal_cohort_chunk_2-calls/"));
        final List<String> chunkDirectoriesFileList = chunkDirectoriesNames.stream()
                .map(s -> TEST_DIR.getAbsolutePath().concat(s)).collect(Collectors.toList());
        final File outputVCF = createTempFile("test", ".vcf");

        final List<String> arguments = new ArrayList<>();
        chunkDirectoriesFileList.forEach(dir -> {
            arguments.add("-" + PostProcessGermlineCNVCalls.CHUNK_PATH_LONG_NAME); arguments.add(dir);
        });
        arguments.add("-" + PostProcessGermlineCNVCalls.SAMPLE_DIRECTORY_LONG_NAME);
        arguments.add(SAMPLE_NAME);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        arguments.add(outputVCF.getAbsolutePath());
        arguments.add("--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE);
        arguments.add("false");

        runCommandLine(arguments);

        IntegrationTestSpec.assertEqualTextFiles(outputVCF, VCF_CORRECT_OUTPUT);
    }
}