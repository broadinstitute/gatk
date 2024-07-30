package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
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

    @DataProvider(name = "getTestFiles")
    public Object[][] getTestFiles() {
        return new Object[][] {
                { " -L chr1:186475", "expected_warning_as_vqslod.txt" },
                { " -L chr1:186475 --ignore-attribute AS_VQSLOD", "empty_file.txt" },
                { " -L chr1:187471 --ignore-attribute AS_VQSLOD", "expected_warning_filter.txt"},
                { " -L chr1:186475-945670 " + DEFAULT_WARP_SETTINGS, "empty_file.txt"}
        };
    }

    @Test(dataProvider = "getTestFiles")
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
        testSpec.executeTest("testDefaultsInWarp", this);
    }

}
