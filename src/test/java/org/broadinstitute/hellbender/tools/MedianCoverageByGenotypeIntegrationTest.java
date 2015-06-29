package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class MedianCoverageByGenotypeIntegrationTest extends CommandLineProgramTest {

    @Test(dataProvider = "filenames")
    public void testBasic(String fileIn) throws Exception {
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        final String[] args = new String[]{
                "--variant" ,  ORIG_FILE.getAbsolutePath(),
                "--output", "testBasic_MedianCoverageByGenotypeIntegrationTest.txt"

        };

        final String res = (String)this.runCommandLine(args);
        System.out.print(res);

        Assert.assertEquals(res, 5);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"CEUTrio-snps.vcf"},
        };
    }
}
