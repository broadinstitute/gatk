package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Testing the test infrastructure
 */
public final class TestSpecUnitTest extends BaseTest {

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

}
