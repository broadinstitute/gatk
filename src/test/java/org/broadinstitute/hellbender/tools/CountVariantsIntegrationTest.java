package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.CountVariants;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountVariantsIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return CountVariants.class.getSimpleName();
    }

    @Test(dataProvider = "filenames")
    public void testCountVariants(final File fileIn, final String moreArgs, final long expected) throws Exception {
        ArgumentsBuilder b= new ArgumentsBuilder();
        b.addArgument("variant", fileIn.getAbsolutePath());
        b.add(moreArgs);
        final Object res = this.runCommandLine(b.getArgsArray());
        Assert.assertEquals(res, expected);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new Object[][]{
                {new File(getTestDataDir(), "count_variants.vcf"), "", 26L},
                {new File(getTestDataDir(), "count_variants.vcf"), "-L 3", 4L},
                {new File(getTestDataDir(), "count_variants.blockgz.gz"), "", 26L},
                {new File(dbsnp_138_b37_1_65M_vcf), "", 1375319L},
        };
    }
}
