package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountVariantsSparkIntegrationTest extends CommandLineProgramTest {

    public static final File COUNT_VARIANTS_VCF = new File(getTestDataDir(), "count_variants.vcf");

    @Override
    public String getTestedClassName() {
        return CountVariantsSpark.class.getSimpleName();
    }

    @Test(dataProvider = "filenames")
    public void test(final File fileIn, final long expected) throws Exception {
        final File outputTxt = createTempFile("count_variants", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(fileIn);
        args.addOutput(outputTxt);
        this.runCommandLine(args.getArgsArray());

        final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile());
        Assert.assertEquals((int)Integer.valueOf(readIn), expected);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new Object[][]{
                {COUNT_VARIANTS_VCF, 26L},
//                {new File(getTestDataDir(), "count_variants.blockgz.gz"), 26L}, //disabled because of https://github.com/HadoopGenomics/Hadoop-BAM/issues/68
                {new File(dbsnp_138_b37_1_65M_vcf), 1375319L},
        };
    }

    @Test
    public void testNoNPRWhenOutputIsUnspecified(){
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(COUNT_VARIANTS_VCF);
        this.runCommandLine(args.getArgsArray());
    }
}