package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.CountVariants;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountVariantsIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return CountVariants.class.getSimpleName();
    }

    @Test(dataProvider = "countVariantsVCFInputs")
    public void testCountVariants(final File fileIn, final String moreArgs, final long expectedCount) throws Exception {
        final ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.addVCF(fileIn);
        ab.add(moreArgs);
        final Object res = runCommandLine(ab.getArgsArray());
        Assert.assertEquals(res, expectedCount);
    }

    @DataProvider(name="countVariantsVCFInputs")
    public Object[][] countVariantsVCFInputs() {
        return new Object[][]{
                {new File(getTestDataDir(), "count_variants.vcf"), "", 26L},
                {new File(getTestDataDir(), "count_variants.blockgz.gz"), "", 26L},
                {new File(getTestDataDir(), "count_variants_withSequenceDict.vcf"), "", 26L},
                {new File(getTestDataDir(), "count_variants_withSequenceDict.vcf"), "-L 1", 14L},
                {new File(dbsnp_138_b37_1_65M_vcf), "", 1375319L},
        };
    }
}
