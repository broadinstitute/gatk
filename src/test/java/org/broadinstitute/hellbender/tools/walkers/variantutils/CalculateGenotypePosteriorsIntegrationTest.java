package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Collections;

public final class CalculateGenotypePosteriorsIntegrationTest extends CommandLineProgramTest {

    private final String largeDir = largeFileTestDir;
    private final String dir = getToolTestDataDir();

    private String CEUtrioFamilyFile = dir + "CEUtrio.ped";
    private String threeMemberNonTrioFamilyFile = dir + "threeMemberNonTrio.ped";

    private String CEUtrioTest = dir + "CEUtrioTest_chr1.vcf";
    private String CEUtrioPopPriorsTest = dir + "CEUtrioPopPriorsTest_chr1.vcf";
    private String CEUtrioMixedPloidyTest = dir + "CEUtrioMixedPloidy.vcf";
    private String getThreeMemberNonTrioTest = dir + "threeMemberNonTrioTest_chr1.vcf";

    @Test
    public void testDefaultsWithPanel() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -O %s" +
                        " -R " + b37_reference_20_21 +    //NOTE: we need a reference for -L
                        " -L 20:10,000,000-10,010,000" +
                        " -" +CalculateGenotypePosteriors.SUPPORTING_CALLSETS_SHORT_NAME + " " + largeDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false" +
                        " -V " + dir + "NA12878.Jan2013.haplotypeCaller.subset.indels.vcf",
                Collections.singletonList(dir + "expectedCGP_testDefaultsWithPanel.vcf")
        );
        spec.executeTest("testDefaultsWithPanel", this);
    }

    @Test
    public void testNumRefWithPanel() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -O %s" +
                        " -R " + b37_reference_20_21 +    //NOTE: we need a reference for -L
                        " -L 20:10,000,000-10,010,000" +
                        " -" +CalculateGenotypePosteriors.SUPPORTING_CALLSETS_SHORT_NAME + " " + largeDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false" +
                        " -V " + dir + "NA12878.Jan2013.haplotypeCaller.subset.indels.vcf" +
                        " --" + CalculateGenotypePosteriors.NUM_REF_SAMPLES_LONG_NAME + " 2500",
                Collections.singletonList(dir + "expectedCGP_testNumRefWithPanel.vcf")
        );
        spec.executeTest("testDefaultsWithPanel", this);
    }

    @Test
    //use the first 20 variants to save time; they have a nice range of AC from 4 to over 4000
    public void testUsingDiscoveredAF() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -O %s" +
                        " -R " + b37_reference_20_21 +    //NOTE: we need a reference for -L
                        " -L 20:10,000,000-10,001,432" +
                        " -V " + largeDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(largeDir + "CalculateGenotypePosteriors/expectedCGP_testUsingDiscoveredAF.vcf")
        );
        spec.executeTest("testUsingDiscoveredAF", this);
    }

    @Test
    //this test ignores discovered AC and has no external priors, so it should only apply the PP tag with values the same as PLs
    //only test the first 20 variants to save time
    public void testMissingPriors() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        "--discovered-allele-count-priors-off" +
                        " -O %s" +
                        " -R " + b37_reference_20_21 +    //NOTE: we need a reference for -L
                        " -L 20:10,000,000-10,001,432" +
                        " -V " + largeDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(largeDir + "CalculateGenotypePosteriors/expectedCGP_testMissingPriors.vcf")
        );
        spec.executeTest("testMissingPriors", this);
    }

    @Test
    public void testInputINDELs() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        "--discovered-allele-count-priors-off" +
                        " -O %s" +
                        " -R " + b37_reference_20_21 +    //NOTE: we need a reference for -L
                        " -L 20:10,000,000-10,100,000" +
                        " -V " + dir + "NA12878.Jan2013.haplotypeCaller.subset.indels.vcf" +
                        " --supporting-callsets " + largeDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(dir + "expectedCGP_testInputINDELs.vcf")
        );
        spec.executeTest("testInputINDELs", this);
    }

    @Test
    public void testFamilyPriors() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        "--discovered-allele-count-priors-off" +
                        " -O %s" +
                        " -ped " + CEUtrioFamilyFile +
                        " -V " + CEUtrioTest +
                        " --supporting-callsets " + CEUtrioPopPriorsTest +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(dir + "expectedCGP_testFamilyPriors_chr1.vcf")
        );
        spec.executeTest("testFamilyPriors", this);
    }

    @Test
    public void testSingleParentFamily() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -O %s" +
                        " -ped " + threeMemberNonTrioFamilyFile +
                        " -V " + getThreeMemberNonTrioTest +
                        " --skip-population-priors" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(dir + "expectedCGP_testSingleParentFamily_chr1.vcf")
        );
        spec.executeTest("testFamilyPriors", this);
    }

    @Test
    public void testFamilyPriorsMixedPloidy() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -O %s" +
                        " -ped " + CEUtrioFamilyFile +
                        " -V " + CEUtrioMixedPloidyTest,
                1,
                UserException.class);
        spec.executeTest("testFamilyPriorsMixedPloidy", this);
    }


    @Test //test for https://github.com/broadinstitute/gatk/issues/4346
    public void testNoDuplicatesForAdjacentInputs() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -O %s" +
                        " -V " + getTestFile("overlappingVariants.vcf") +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE +" false",
                Collections.singletonList(getTestFile("overlappingVariants.expected.no_duplicates.vcf").toString()));
        spec.executeTest("testNoDuplicatesForOverlappingVariants", this);
    }

    @Test //test for https://github.com/broadinstitute/gatk/issues/4346
    public void gzipOutputIsGzipped() throws IOException {
        final File out = createTempFile("out", ".vcf.gz");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(out)
            .addVCF(getTestFile("overlappingVariants.vcf"));

        runCommandLine(args);

        try( final InputStream in = new BufferedInputStream(new FileInputStream(out))) {
            Assert.assertTrue(BlockCompressedInputStream.isValidFile(in));
        }
    }

    @DataProvider
    public Object[][] getBadInputs(){
        return new Object[][]{
                {"noGenotypes.vcf"},
                {"badMLEAC_count_not_A.vcf"},
                {"badMLEAC_type_not_int.vcf"}
        };
    }

    @Test(dataProvider = "getBadInputs", expectedExceptions = UserException.BadInput.class)
    public void testBadInputFilesAreRejectedWithReasonableError(String badFile) throws IOException {
        final File out = createTempFile("out", ".vcf.gz");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(out)
            .addVCF(getTestFile(badFile));

        runCommandLine(args);
    }
}
