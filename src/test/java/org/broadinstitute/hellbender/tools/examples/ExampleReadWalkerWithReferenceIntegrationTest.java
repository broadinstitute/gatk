package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public final class ExampleReadWalkerWithReferenceIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = exampleTestDir;

    @DataProvider
    public Object[][] getReferences(){
        return new Object[][]{
                {hg19MiniReference},
                {hg19MiniReference + ".gz"}
        };
    }


    @Test(dataProvider = "getReferences")
    public void testExampleReadWalkerWithReference(String reference) throws IOException {
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + reference +
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expected_ExampleReadWalkerWithReferenceIntegrationTest_output.txt")
        );
        testSpec.executeTest("testExampleReadWalkerWithReference", this);
    }



}
