package org.broadinstitute.hellbender.testutils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.PrintReadsIntegrationTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Testing the test infrastructure
 */
public final class IntegrationTestSpecUnitTest extends GATKBaseTest {

    @Test(expectedExceptions = AssertionError.class)
    public void compareTextFiles() throws IOException {
        IntegrationTestSpec.assertEqualTextFiles(new File("AUTHORS"), new File("README.md"));
    }

    @Test(expectedExceptions = AssertionError.class)
    public void compareTextFiles_differentEvenIgnoringComments() throws IOException {
        IntegrationTestSpec.assertEqualTextFiles(new File("AUTHORS"), new File("README.md"), "#");
    }

    @Test
    public void compareEqualTextFiles() throws IOException {
        IntegrationTestSpec.assertEqualTextFiles(new File("AUTHORS"), new File("AUTHORS"));
    }

    @DataProvider(name = "differentFilesButSameContent")
    public Object[][] differentFilesButSameContent(){
        return new Object[][]{
                {"fileWithComments.txt", "fileWithNoComments.txt"},
                {"fileWithComments.txt", "fileWithDifferentComments.txt"},
                {"fileWithNoComments.txt", "fileWithComments.txt"},
                {"fileWithNoComments.txt", "fileWithDifferentComments.txt"},
        };
    }

    @Test(expectedExceptions = AssertionError.class, dataProvider = "differentFilesButSameContent")
    public void compareFilesWithAndWithoutComments(final String f1, final String f2) throws IOException {
        final String dir = publicTestDir + "org/broadinstitute/hellbender/utils/testing/";
        IntegrationTestSpec.assertEqualTextFiles(new File(dir + f1), new File(dir + f2));
    }

    @Test(dataProvider = "differentFilesButSameContent")
    public void compareFilesWithAndWithoutComments_ignoreComments(final String f1, final String f2) throws IOException {
        final String dir = publicTestDir + "org/broadinstitute/hellbender/utils/testing/";
        IntegrationTestSpec.assertEqualTextFiles(new File(dir + f1), new File(dir + f2), "#");
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSpec_misspecified() throws IOException {
        final File getTestDataDir = new File(toolsTestDir);
        final File samWithOneMalformedRead = new File(getTestDataDir, "print_reads_one_malformed_read.sam");
        final File outBam = createTempFile("print_reads_testReadFiltering", ".bam");

        //This test spec is incorrect because has 1 expected file and no placeholders for output files
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --" + StandardArgumentDefinitions.INPUT_LONG_NAME + " " + samWithOneMalformedRead.getCanonicalPath() +
                " --" + StandardArgumentDefinitions.OUTPUT_LONG_NAME + " " + outBam.getCanonicalPath(),
                Arrays.asList(new File(getTestDataDir, "expected.print_reads_one_malformed_read.bam").getCanonicalPath())
        );
        spec.executeTest("testReadFiltering_asIntegrationTestSpec", new PrintReadsIntegrationTest());
    }

}
