package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.CountVariants;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;

public final class CountVariantsIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return CountVariants.class.getSimpleName();
    }

    @Test(dataProvider = "countVariantsVCFInputs")
    public void testCountVariants(final File fileIn, final String moreArgs, final long expectedCount) throws Exception {
        final ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.addVCF(fileIn);
        ab.addRaw(moreArgs);
        final Object res = runCommandLine(ab);
        Assert.assertEquals(res, expectedCount);
    }

    @Test(dataProvider = "countVariantsVCFInputs")
    public void testCountVariantsWithOutputFile(final File fileIn, final String moreArgs, final long expectedCount) throws Exception {
        final File output = createTempFile("testCountVariantsWithOutputFile", ".txt");

        final ArgumentsBuilder ab = new ArgumentsBuilder();
        ab.addVCF(fileIn);
        ab.addOutput(output);
        ab.addRaw(moreArgs);

        final Object res = runCommandLine(ab);
        Assert.assertEquals(res, expectedCount);

        Assert.assertEquals(Files.readAllBytes(output.toPath()), Long.toString(expectedCount).getBytes());
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
