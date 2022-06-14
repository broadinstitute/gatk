package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

@Test(groups = {"variantcalling"})
public class RampedHaplotypeCallerIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final String TEST_FILES_DIR = toolsTestDir + "haplotypecaller/";

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasicNoRamps() throws Exception {
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "ramps/test_noramps.expected.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = buildCommonArguments(input, outputPath);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @DataProvider(name = "offrampTestData")
    public Object[][] getOfframpTestData() {

        final Object[][]        testData = {
                { "PRE_FILTER_OFF", "ramps/test_pre_filter_offramp.expected.zip" },
                { "PRE_ASSEMBLER_OFF", "ramps/test_pre_assembler_offramp.expected.zip" },
                { "POST_ASSEMBLER_OFF", "ramps/test_post_assembler_offramp.expected.zip" }
        };

        return testData;
    }

    // THis test has been disabled due to issues with the ZipFilesystem comparisons operations.
    @Test(dataProvider = "offrampTestData", enabled = false)
    public void testOfframp(final String type, final String expectedFilename) throws Exception {
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");
        final File offramp = createTempFile("offramp", ".zip");

        final File offrampExpected = new File(TEST_FILES_DIR, expectedFilename);
        final String outputPath = output.getAbsolutePath();
        final String offrampPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? offrampExpected.getAbsolutePath() : offramp.getAbsolutePath();

        final ArgumentsBuilder args = buildCommonArguments(input, outputPath);
        args.add("off-ramp-type", type);
        args.add("off-ramp-file", offrampPath);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualZipFiles(offramp, offrampExpected, IntegrationTestSpec.EqualZipFilesAssist_AllText);
        }
    }

    @DataProvider(name = "onrampTestData")
    public Object[][] getOnrampTestData() {

        final Object[][]        testData = {
                { "POST_ASSEMBLER_ON", "ramps/test_post_assembler_offramp.expected.zip", "ramps/test_post_assembler_output.expected.vcf" }
        };

        return testData;
    }

    @Test(dataProvider = "onrampTestData")
    public void testOnramp(final String type, final String rampFilename, final String expectedFilename) throws Exception {
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, expectedFilename);
        final File onrampFile = new File(TEST_FILES_DIR, rampFilename);
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = buildCommonArguments(input, outputPath);
        args.add("on-ramp-type", type);
        args.add("on-ramp-file", onrampFile);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    private ArgumentsBuilder buildCommonArguments(File input, String outputPath) {

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBased")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .add("flow-likelihood-optimized-comp", true)
                .add("ERC", "GVCF")
                .add("flow-filter-lone-alleles", true)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        return args;
    }
}
