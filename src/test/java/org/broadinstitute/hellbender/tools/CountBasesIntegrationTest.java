package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class CountBasesIntegrationTest extends CommandLineProgramTest {

    @Test(dataProvider = "filenames")
    public void testCountBases(String fileIn) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final String[] args = new String[]{
            "--input" ,  ORIG_BAM.getAbsolutePath(),

        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, 808l);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"count_bases.sam"},
                {"count_bases.bam"},
        };
    }
}
