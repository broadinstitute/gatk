package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class VCFComparatorIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/walkers/VCFComparator/";
    private static final String TEST_OUTPUT_DIRECTORY = TEST_DATA_DIRECTORY + "/expected/";
    private static final String DEFAULT_WARP_SETTINGS = "--ignore-attribute VQSLOD --ignore-attribute AS_VQSLOD --ignore-filters --ignore-attribute culprit  --ignore-attribute AS_culprit --ignore-attribute AS_FilterStatus --ignore-attribute ExcessHet --ignore-attribute AS_SOR --ignore-attribute AS_FS --ignore-attribute AS_BaseQRankSum --ignore-attribute AS_ReadPosRankSum --ignore-attribute AS_MQRankSum";

    @DataProvider(name = "getTestVcfs")
    public Object[][] getTestVcfs() {
        return new Object[][] {
                { " -L chr1:186475", "expected_warning_as_vqslod.txt" },
                { " -L chr1:186475 --ignore-attribute AS_VQSLOD", "empty_file.txt" },
                { " -L chr1:187471 --ignore-attribute AS_VQSLOD", "expected_warning_filter.txt" },
                { " -L chr1:186475-945669 " + DEFAULT_WARP_SETTINGS, "empty_file.txt" },
                { " -L chr1:945670 " + DEFAULT_WARP_SETTINGS, "qual_diff_warning.txt"}, // different QUAL values
                { " -L chr1:945670 --qual-change-allowed 0.1 --ignore-attribute AS_VQSLOD", "empty_file.txt"},
                { " -L chr1:186475 --mute-acceptable-diffs", "empty_file.txt" } // low quality site is muted even though the AS_VQSLOD is different
        };
    }

    @Test(dataProvider = "getTestVcfs")
    public void testAnnotationDifferences(String args, String expectedWarnings) throws IOException {
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                        " -R " + hg38Reference +
                        " -V:actual " + TEST_DATA_DIRECTORY + "actual.vcf" +
                        " -V:expected " + TEST_DATA_DIRECTORY + "expected.vcf" +
                        " --output-warnings %s" +
                        " --warn-on-errors" +
                        args,
                Arrays.asList(TEST_OUTPUT_DIRECTORY + expectedWarnings)
        );
        testSpec.executeTest("testVcfs", this);
    }

    @Test
    public void testExpectedFailure() {
        final File actual = new File(TEST_DATA_DIRECTORY + "actual.vcf");
        final File expected = new File(TEST_DATA_DIRECTORY + "expected.vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(hg38Reference);
        args.add("V:actual", actual);
        args.add("V:expected", expected);
        args.addInterval("chr1:186475");
        args.add("warn-on-errors", true);

        //warn on errors alone shouldn't fail
        runCommandLine(args);

        //adding finish-before-failing should fail
        args.add("finish-before-failing", true);
        Assert.assertThrows(UserException.class, () -> runCommandLine(args));
    }

    @DataProvider(name = "getTestGvcfs")
    public Object[][] getTestGvcfs() {
        return new Object[][] {
                { " -L chr1:864084-864610", "empty_file.txt" }, //matching ref blocks
                { " -L chr1:54682-347969", "ref_block_warning.txt"}, // non-matching ref block
                { " -L chr1:792417", "tree_score_warning.txt"}, // variant site
                { " -L chr1:792417 --ignore-non-ref-data --" +
                        ReblockGVCF.ANNOTATIONS_TO_KEEP_LONG_NAME + " TREE_SCORE", "tree_score_warning.txt"}, // when non-ref data is dropped non-GATK annotations can be dropped
                { " -L chr1:792417 --ignore-attribute TREE_SCORE", "empty_file.txt"}
        };
    }

    @Test(dataProvider = "getTestGvcfs")
    public void testGvcfs(String args, String expectedWarnings) throws IOException {
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + hg38Reference +
                        " -V:actual " + TEST_DATA_DIRECTORY + "actual.NA12878.rb.g.vcf" +
                        " -V:expected " + TEST_DATA_DIRECTORY + "expected.NA12878.rb.g.vcf" +
                        " --output-warnings %s" +
                        " --warn-on-errors" +
                        args,
                Arrays.asList(TEST_OUTPUT_DIRECTORY + expectedWarnings));
        testSpec.executeTest("testGvcfs", this);
    }

    @Test
    public void testMixedPloidy() throws IOException {
        final IntegrationTestSpec expectDiploid = new IntegrationTestSpec(
                " -R " + hg38Reference +
                        " -L chrX:66780645" +
                        " -V:actual " + TEST_DATA_DIRECTORY + "haploid.rb.g.vcf" +
                        " -V:expected " + TEST_DATA_DIRECTORY + "diploid.rb.g.vcf" +
                        " --output-warnings %s" +
                        " --warn-on-errors",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "haploid_warning.txt"));
        expectDiploid.executeTest("testMixedPloidy", this);

        // check the other way as far as which ploidy is expected vs actual
        final IntegrationTestSpec expectHaploid = new IntegrationTestSpec(
                " -R " + hg38Reference +
                        " -L chrX:66780645" +
                        " -V:actual " + TEST_DATA_DIRECTORY + "diploid.rb.g.vcf" +
                        " -V:expected " + TEST_DATA_DIRECTORY + "haploid.rb.g.vcf" +
                        " --output-warnings %s" +
                        " --warn-on-errors",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "haploid_warning2.txt"));
        expectHaploid.executeTest("testMixedPloidy", this);
    }

}
