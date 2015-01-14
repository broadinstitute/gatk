package org.broadinstitute.hellbender.tools.picard;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;

/**
 * ValidateSamFile is a thin wrapper around {@link htsjdk.samtools.SamFileValidator} which is thoroughly tested in HTSJDK.
 * Thus we don't have much to test here.
 */
public class ValidateSamFileTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/ValidateSamFile");

    @Override
    public String getCommandLineProgramName() {
        return ValidateSamFile.class.getSimpleName();
    }

    @DataProvider(name="testingData")
    public Object[][] testingData() {
        return new Object[][]{
            {"valid.sam", "SUMMARY", true},
            {"valid.sam", "VERBOSE", true},
            {"invalid_coord_sort_order.sam", "SUMMARY", false},
        };
    }

    @Test(dataProvider="testingData")
    public void testNoOutputFile(String input, String mode, boolean expectedValidity) throws Exception {
        final File testFile = new File(TEST_DATA_DIR, input);
        final String[] args = new String[] {
            "--INPUT", testFile.getPath(),
            "--MODE", mode
        };
        Assert.assertEquals(runCommandLine(args), expectedValidity);
    }
}
