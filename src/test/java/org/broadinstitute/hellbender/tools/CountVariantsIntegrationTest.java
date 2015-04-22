package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountVariantsIntegrationTest extends CommandLineProgramTest {

    @Test(dataProvider = "filenames")
    public void testCountBases(String fileIn) throws Exception {
        final File ORIG_FILE = new File(getTestDataDir(), fileIn);
        final String[] args = new String[]{
            "--variant" ,  ORIG_FILE.getAbsolutePath(),

        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, 26l);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"count_variants.vcf"},
        };
    }
}
