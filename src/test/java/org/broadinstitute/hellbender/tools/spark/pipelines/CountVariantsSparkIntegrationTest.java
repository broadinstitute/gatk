package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

public final class CountVariantsSparkIntegrationTest extends CommandLineProgramTest {

    public static final File COUNT_VARIANTS_VCF = new File(getTestDataDir(), "count_variants.vcf");

    @Override
    public String getTestedClassName() {
        return CountVariantsSpark.class.getSimpleName();
    }

    @Test(dataProvider = "filenames", groups = "spark")
    public void test(final File fileIn, final long expected) throws Exception {
        final File outputTxt = createTempFile("count_variants", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(fileIn);
        args.addOutput(outputTxt);
        this.runCommandLine(args.getArgsArray());

        final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile(), StandardCharsets.UTF_8);
        Assert.assertEquals((int)Integer.valueOf(readIn), expected);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new Object[][]{
                {COUNT_VARIANTS_VCF, 26L},
                {new File(getTestDataDir(), "count_variants.blockgz.gz"), 26L},
                {new File(dbsnp_138_b37_1_65M_vcf), 1375319L},
        };
    }

    @DataProvider(name="intervals")
    public Object[][] intervals(){
        File vcf = new File(largeFileTestDir, "dbsnp_138.b37.20.21.vcf");
        File vcf_gz = new File(largeFileTestDir, "dbsnp_138.b37.20.21.vcf.blockgz.gz");
        return new Object[][]{
                new Object[]{vcf, "", 9594L}, // no intervals specified
                new Object[]{vcf, "-L 20", 5657L},
                new Object[]{vcf, "-L 20:10200000-11000000", 933L},
                new Object[]{vcf, "-L 21", 3937L},
                new Object[]{vcf, "-L 20 -L 21", 9594L},
                new Object[]{vcf, "-XL 20", 3937L},
                new Object[]{vcf_gz, "", 9594L}, // no intervals specified
                new Object[]{vcf_gz, "-L 20", 5657L},
                new Object[]{vcf_gz, "-L 20:10200000-11000000", 933L},
                new Object[]{vcf_gz, "-L 21", 3937L},
                new Object[]{vcf_gz, "-L 20 -L 21", 9594L},
                new Object[]{vcf_gz, "-XL 20", 3937L},
        };
    }

    @Test(dataProvider = "intervals", groups = "spark")
    public void testCountVariantsWithIntervals(final File fileIn, final String intervalArgs, final long expected) throws Exception {
        final File outputTxt = createTempFile("count_variants", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(fileIn);
        args.add(intervalArgs);
        args.addReference(new File(largeFileTestDir, "human_g1k_v37.20.21.fasta"));
        args.addOutput(outputTxt);
        this.runCommandLine(args.getArgsArray());

        final ByteArrayOutputStream baosErr = new ByteArrayOutputStream();
        final PrintStream err = System.err;
        try {
            System.setErr(new PrintStream(baosErr));

            this.runCommandLine(args.getArgsArray());

            final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile(), StandardCharsets.UTF_8);
            Assert.assertEquals((int)Integer.valueOf(readIn), expected);
            String errString = baosErr.toString();
            Assert.assertFalse(errString.contains("Warning: using GzipCodec, which is not splittable,"), errString);
        } finally {
            System.setErr(err); //put this back in
        }
    }

    @Test(groups = "spark")
    public void testNoNPRWhenOutputIsUnspecified(){
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(COUNT_VARIANTS_VCF);
        this.runCommandLine(args.getArgsArray());
    }
}
