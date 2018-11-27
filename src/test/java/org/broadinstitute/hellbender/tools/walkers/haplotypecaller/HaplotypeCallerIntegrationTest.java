package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SamFiles;
import htsjdk.tribble.Tribble;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AssemblyRegionWalker;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

@Test(groups = {"variantcalling"})
public class HaplotypeCallerIntegrationTest extends AbstractHaplotypeCallerIntegrationTest {

    /*
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }


    /*
     * Test that in GVCF mode we're consistent with past GATK4 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testGVCFModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults", ".g.vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testGVCFMode.gatk4.g.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that in GVCF mode we're consistent with past GATK4 results using AS_ annotations
     *
     * Updated on 09/01/17 to account for changes to AS_RankSum annotations the annotations were checked against GATK3
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testGVCFModeIsConsistentWithPastResults_AlleleSpecificAnnotations(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults_AlleleSpecificAnnotations", ".g.vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testGVCFMode.gatk4.alleleSpecific.g.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @Test
    public void testGVCFModeGenotypePosteriors() throws Exception {
        Utils.resetRandomGenerator();

        final String inputFileName = NA12878_20_21_WGS_bam;
        final String referenceFileName =b37_reference_20_21;

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults", ".g.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ERC", "GVCF",
                "--" + GenotypeCalculationArgumentCollection.SUPPORTING_CALLSET_LONG_NAME,
                    largeFileTestDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf",
                "--" + GenotypeCalculationArgumentCollection.NUM_REF_SAMPLES_LONG_NAME, "2500",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        Pair<VCFHeader, List<VariantContext>> results = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        for (final VariantContext vc : results.getRight()) {
            final Genotype g = vc.getGenotype(0);
            if (g.hasDP() && g.getDP() > 0 && g.hasGQ() && g.getGQ() > 0) {
                Assert.assertTrue(g.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
            }
            if (isGVCFReferenceBlock(vc) ) {
                Assert.assertTrue(!vc.hasAttribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
            }
            else if (!vc.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE)){      //there are some variants that don't have non-symbolic alts
                Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
            }
        }
    }

    @Test
    public void testGenotypeGivenAllelesMode() throws IOException {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGenotypeGivenAllelesMode", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testGenotypeGivenAllelesMode.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                "--genotyping-mode", "GENOTYPE_GIVEN_ALLELES",
                "--alleles", new File(TEST_FILES_DIR, "testGenotypeGivenAllelesMode_givenAlleles.vcf").getAbsolutePath()
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testGenotypeGivenAllelesModeNotAllowedInGVCFMode() throws IOException {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGenotypeGivenAllelesModeNotAllowedInGVCFMode", ".g.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                "--genotyping-mode", "GENOTYPE_GIVEN_ALLELES",
                "--alleles", new File(TEST_FILES_DIR, "testGenotypeGivenAllelesMode_givenAlleles.vcf").getAbsolutePath(),
                "-ERC", "GVCF"
        };

        // Should throw, since -ERC GVCF is incompatible with GENOTYPE_GIVEN_ALLELES mode
        runCommandLine(args);
    }

    @Test
    public void testBamoutProducesReasonablySizedOutput() {
        final Path bamOutput = createTempFile("testBamoutProducesReasonablySizedOutput", ".bam").toPath();
        innerTestBamoutProducesReasonablySizedOutput(bamOutput);
    }

    @Test(groups={"bucket"})
    public void testBamoutOnGcs() {
        final Path bamOutput = BucketUtils.getPathOnGcs(BucketUtils.getTempFilePath(
            getGCPTestStaging() + "testBamoutProducesReasonablySizedOutput", ".bam"));
        innerTestBamoutProducesReasonablySizedOutput(bamOutput);
    }

    private void innerTestBamoutProducesReasonablySizedOutput(Path bamOutput) {
        Utils.resetRandomGenerator();

        // We will test that when running with -bamout over the testInterval, we produce
        // a bam with a number of reads that is within 10% of what GATK3.5 produces with
        // -bamout over the same interval. This is just to test that we produce a reasonably-sized
        // bam for the region, not to validate the haplotypes, etc. We don't want
        // this test to fail unless there is a likely problem with -bamout itself (eg., empty
        // or truncated bam).
        final String testInterval = "20:10000000-10010000";
        final int gatk3BamoutNumReads = 5170;

        final File vcfOutput = createTempFile("testBamoutProducesReasonablySizedOutput", ".vcf");

        ArgumentsBuilder argBuilder = new ArgumentsBuilder();

        argBuilder.addInput(new File(NA12878_20_21_WGS_bam));
        argBuilder.addReference(new File(b37_reference_20_21));
        argBuilder.addOutput(new File(vcfOutput.getAbsolutePath()));
        argBuilder.addArgument("L", testInterval);
        argBuilder.addArgument(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_SHORT_NAME, bamOutput.toUri().toString());
        argBuilder.addArgument("pairHMM", "AVX_LOGLESS_CACHING");

        runCommandLine(argBuilder.getArgsArray());

        try ( final ReadsDataSource bamOutReadsSource = new ReadsDataSource(bamOutput) ) {
            int actualBamoutNumReads = 0;
            for ( final GATKRead read : bamOutReadsSource ) {
                ++actualBamoutNumReads;
            }

            final int readCountDifference = Math.abs(actualBamoutNumReads - gatk3BamoutNumReads);
            Assert.assertTrue(((double)readCountDifference / gatk3BamoutNumReads) < 0.10,
                    "-bamout produced a bam with over 10% fewer/more reads than expected");
        }
    }

    @Test
    public void testSitesOnlyMode() {
        Utils.resetRandomGenerator();
        File out = createTempFile("GTStrippedOutput", "vcf");
        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", out.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.SITES_ONLY_LONG_NAME,
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Assert that the genotype field has been stripped from the file
        Pair<VCFHeader, List<VariantContext>> results = VariantContextTestUtils.readEntireVCFIntoMemory(out.getAbsolutePath());

        Assert.assertFalse(results.getLeft().hasGenotypingData());
        for (VariantContext v: results.getRight()) {
            Assert.assertFalse(v.hasGenotypes());
        }
    }

    @Test
    public void testForceActiveOption() throws Exception {
        Utils.resetRandomGenerator();
        final File out = createTempFile("GTStrippedOutput", "vcf");
        final File assemblyRegionOut = createTempFile("assemblyregions", ".igv");
        final File expectedAssemblyRegionOut = new File(TEST_FILES_DIR, "expected.testAssemblyRegionWithForceActiveRegions_assemblyregions.igv");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:1-5000",
                "-O", out.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyRegionWalker.FORCE_ACTIVE_REGIONS_LONG_NAME, "true",
                "--" + HaplotypeCaller.ASSEMBLY_REGION_OUT_LONG_NAME, assemblyRegionOut.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(assemblyRegionOut, expectedAssemblyRegionOut);
    }

    @DataProvider(name="outputFileVariations")
    public Object[][] getOutputFileVariations() {
        return new Object[][]{
                // bamout index, bamout md5, vcf index, vcf md5
                { true, true, true, true },
                { true, false, true, false },
                { false, true, false, true },
                { false, false, false, false },
        };
    }

    @Test(dataProvider = "outputFileVariations")
    public void testOutputFileArgumentVariations(
            final boolean createBamoutIndex,
            final boolean createBamoutMD5,
            final boolean createVCFOutIndex,
            final boolean createVCFOutMD5) throws IOException {
        Utils.resetRandomGenerator();

        // run on small interval to test index/md5 outputs
        final String testInterval = "20:10000000-10001000";

        final File vcfOutput = createTempFile("testOutputFileArgumentVariations", ".vcf");
        final File bamOutput = createTempFile("testOutputFileArgumentVariations", ".bam");

        ArgumentsBuilder argBuilder = new ArgumentsBuilder();

        argBuilder.addInput(new File(NA12878_20_21_WGS_bam));
        argBuilder.addReference(new File(b37_reference_20_21));
        argBuilder.addOutput(new File(vcfOutput.getAbsolutePath()));
        argBuilder.addArgument("L", testInterval);
        argBuilder.addArgument(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_SHORT_NAME, bamOutput.getAbsolutePath());
        argBuilder.addArgument("pairHMM", "AVX_LOGLESS_CACHING");
        argBuilder.addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_LONG_NAME, createBamoutIndex);
        argBuilder.addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_LONG_NAME, createBamoutMD5);
        argBuilder.addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME, createVCFOutIndex);
        argBuilder.addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_MD5_LONG_NAME, createVCFOutMD5);

        runCommandLine(argBuilder.getArgsArray());

        Assert.assertTrue(vcfOutput.exists(), "No VCF output file was created");

        // validate vcfout companion files
        final File vcfOutFileIndex = new File(vcfOutput.getAbsolutePath() + Tribble.STANDARD_INDEX_EXTENSION);
        final File vcfOutFileMD5 = new File(vcfOutput.getAbsolutePath() + ".md5");
        Assert.assertEquals(vcfOutFileIndex.exists(), createVCFOutIndex, "The index file argument was not honored");
        Assert.assertEquals(vcfOutFileMD5.exists(), createVCFOutMD5, "The md5 file argument was not honored");

        // validate bamout companion files
        if (createBamoutIndex) {
            Assert.assertNotNull(SamFiles.findIndex(bamOutput));
        } else {
            Assert.assertNull(SamFiles.findIndex(bamOutput));
        }

        final File expectedBamoutMD5File = new File(bamOutput.getAbsolutePath() + ".md5");
        Assert.assertEquals(expectedBamoutMD5File.exists(), createBamoutMD5);

        // Check the output BAN header contains all of the inout BAM header Program Records (@PG)
        SamAssertionUtils.assertOutBamContainsInBamProgramRecords(new File(NA12878_20_21_WGS_bam), bamOutput);
    }

    @Test
    public void testAssemblyRegionAndActivityProfileOutput() throws Exception {
        final File output = createTempFile("testAssemblyRegionAndActivityProfileOutput", ".vcf");
        final File assemblyRegionOut = createTempFile("testAssemblyRegionAndActivityProfileOutput_assemblyregions", ".igv");
        final File activityProfileOut = createTempFile("testAssemblyRegionAndActivityProfileOutput_activityprofile", ".igv");
        final File expectedAssemblyRegionOut = new File(TEST_FILES_DIR, "expected.testAssemblyRegionAndActivityProfileOutput_assemblyregions.igv");
        final File expectedActivityProfileOut = new File(TEST_FILES_DIR, "expected.testAssemblyRegionAndActivityProfileOutput_activityprofile.igv");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10003000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + HaplotypeCaller.ASSEMBLY_REGION_OUT_LONG_NAME, assemblyRegionOut.getAbsolutePath(),
                "--" + HaplotypeCaller.PROFILE_OUT_LONG_NAME, activityProfileOut.getAbsolutePath()
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(assemblyRegionOut, expectedAssemblyRegionOut);
        IntegrationTestSpec.assertEqualTextFiles(activityProfileOut, expectedActivityProfileOut);
    }

    // test on an artificial bam with several contrived MNPs
    // this test is basically identical to a test in {@ link Mutect2IntegrationTest}
    @Test
    public void testMnps() throws Exception {
        Utils.resetRandomGenerator();
        final File bam = new File(toolsTestDir, "mnp.bam");

        for (final int maxMnpDistance : new int[] {0, 1, 2, 3, 5}) {
            final File outputVcf = createTempFile("unfiltered", ".vcf");

            final List<String> args = Arrays.asList("-I", bam.getAbsolutePath(),
                    "-R", b37_reference_20_21,
                    "-L", "20:10019000-10022000",
                    "-O", outputVcf.getAbsolutePath(),
                    "-" + HaplotypeCallerArgumentCollection.MAX_MNP_DISTANCE_SHORT_NAME, Integer.toString(maxMnpDistance));
            runCommandLine(args);

            checkMnpOutput(maxMnpDistance, outputVcf);
        }
    }

    // this is particular to our particular artificial MNP bam -- we extract a method in order to use it for HaplotypeCaller
    private static void checkMnpOutput(int maxMnpDistance, File outputVcf) {
        // note that for testing HaplotypeCaller GVCF mode we will always have the symbolic <NON REF> allele
        final Map<Integer, List<String>> alleles = VariantContextTestUtils.streamVcf(outputVcf)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc.getAlternateAlleles().stream().filter(a -> !a.isSymbolic()).map(Allele::getBaseString).collect(Collectors.toList())));

        // phased, two bases apart
        if (maxMnpDistance < 2) {
            Assert.assertEquals(alleles.get(10019968), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10019970), Arrays.asList("G"));
        } else {
            Assert.assertEquals(alleles.get(10019968), Arrays.asList("GAG"));
            Assert.assertTrue(!alleles.containsKey(10019970));
        }

        // adjacent and out of phase
        Assert.assertEquals(alleles.get(10020229), Arrays.asList("A"));
        Assert.assertEquals(alleles.get(10020230), Arrays.asList("G"));

        // 4-substitution MNP w/ spacings 2, 3, 4
        if (maxMnpDistance < 2) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020432), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020435), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else if (maxMnpDistance < 3) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAG"));
            Assert.assertEquals(alleles.get(10020435), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else if (maxMnpDistance < 4) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAGTTG"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAGTTGTCTG"));
        }

        // two out of phase DNPs that overlap and have a base in common
        if (maxMnpDistance > 0) {
            Assert.assertEquals(alleles.get(10020680), Arrays.asList("TA"));
            Assert.assertEquals(alleles.get(10020681), Arrays.asList("AT"));
        }
    }

    @DataProvider
    public Object[][] getMaxAlternateAllelesData() {
        return new Object[][] {
                // bam, reference, interval string, max alternate alleles, GVCF mode toggle
                { NA12878_20_21_WGS_bam, b37_reference_20_21, "20:10008000-10010000", 1, false },
                { NA12878_20_21_WGS_bam, b37_reference_20_21, "20:10002000-10011000", 1, true },
                { NA12878_20_21_WGS_bam, b37_reference_20_21, "20:10002000-10011000", 2, true }
        };
    }

    /*
     * Test for the --max-alternate-alleles argument
     */
    @Test(dataProvider = "getMaxAlternateAllelesData")
    public void testMaxAlternateAlleles(final String bam, final String reference, final String intervalString,
                                        final int maxAlternateAlleles, final boolean gvcfMode) {
        final File outputNoMaxAlternateAlleles = createTempFile("testMaxAlternateAllelesNoMaxAlternateAlleles", (gvcfMode ? ".g.vcf" : ".vcf"));
        final File outputWithMaxAlternateAlleles = createTempFile("testMaxAlternateAllelesWithMaxAlternateAlleles", (gvcfMode ? ".g.vcf" : ".vcf"));

        // Run both with and without --max-alternate-alleles over our interval, so that we can
        // prove that the argument is working as intended.
        final String[] argsNoMaxAlternateAlleles = {
                "-I", bam,
                "-R", reference,
                "-L", intervalString,
                "-O", outputNoMaxAlternateAlleles.getAbsolutePath(),
                "-ERC", (gvcfMode ? "GVCF" : "NONE")
        };
        runCommandLine(argsNoMaxAlternateAlleles);

        final String[] argsWithMaxAlternateAlleles = {
                "-I", bam,
                "-R", reference,
                "-L", intervalString,
                "-O", outputWithMaxAlternateAlleles.getAbsolutePath(),
                "--max-alternate-alleles", Integer.toString(maxAlternateAlleles),
                "-ERC", (gvcfMode ? "GVCF" : "NONE")
        };
        runCommandLine(argsWithMaxAlternateAlleles);

        final List<VariantContext> callsNoMaxAlternateAlleles = VariantContextTestUtils.readEntireVCFIntoMemory(outputNoMaxAlternateAlleles.getAbsolutePath()).getRight();
        final List<VariantContext> callsWithMaxAlternateAlleles = VariantContextTestUtils.readEntireVCFIntoMemory(outputWithMaxAlternateAlleles.getAbsolutePath()).getRight();

        // First, find all calls in the VCF produced WITHOUT --max-alternate-alleles that have
        // more than maxAlternateAlleles alt alleles, excluding NON_REF. For each call, calculate
        // and store the expected list of subsetted alleles:
        final Map<SimpleInterval, List<Allele>> expectedSubsettedAllelesByLocus = new HashMap<>();
        for ( final VariantContext vc : callsNoMaxAlternateAlleles ) {
            if ( getNumAltAllelesExcludingNonRef(vc) > maxAlternateAlleles ) {
                final List<Allele> mostLikelyAlleles = AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, HomoSapiensConstants.DEFAULT_PLOIDY, maxAlternateAlleles);
                expectedSubsettedAllelesByLocus.put(new SimpleInterval(vc), mostLikelyAlleles);
            }
        }

        // Then assert that we saw at least one call with more than maxAlternateAlleles alt alleles
        // when running without --max-alternate-alleles (otherwise, the tests below won't be meaningful):
        Assert.assertTrue(! expectedSubsettedAllelesByLocus.isEmpty(),
                "Without --max-alternate-alleles, there should be at least one call in the output with more than " + maxAlternateAlleles +
                        " alt alleles in order for this test to be meaningful");

        // Now assert that in the VCF produced WITH --max-alternate-alleles, there are no calls with
        // more than maxAlternateAlleles alt alleles, excluding NON_REF. Also check each call that would
        // have had more than maxAlternateAlleles alleles against the expected list of subsetted alleles,
        // to ensure that we selected the most likely alleles:
        for ( final VariantContext vc : callsWithMaxAlternateAlleles ) {

            // No call should have more than the configured number of alt alleles (excluding NON_REF)
            Assert.assertTrue(getNumAltAllelesExcludingNonRef(vc) <= maxAlternateAlleles,
                    "Number of alt alleles exceeds --max-alternate-alleles " + maxAlternateAlleles + " for VariantContext: " + vc);

            // If there's an entry for this locus in our table of expected alleles post-subsetting, assert
            // that we selected the right alleles during subsetting.
            List<Allele> alleleSubsettingExpectedResult = expectedSubsettedAllelesByLocus.get(new SimpleInterval(vc));
            if ( alleleSubsettingExpectedResult != null ) {

                // CollectionUtils.isEqualCollection() will compare the lists of Alleles without
                // regard to ordering
                Assert.assertTrue(CollectionUtils.isEqualCollection(vc.getAlleles(), alleleSubsettingExpectedResult),
                        "For call " + vc + " expected alleles after subsetting were: " + alleleSubsettingExpectedResult +
                                 " but instead found alleles: " + vc.getAlleles());
            }

            // For completeness sake, also check the genotypes to ensure that no genotypes reference
            // an allele not present in the VC:
            for ( final Genotype genotype : vc.getGenotypes() ) {
                if ( genotype.isAvailable() ) {
                    for ( final Allele genotypeAllele : genotype.getAlleles() ) {
                        if ( genotypeAllele.isCalled() ) {
                            Assert.assertTrue(vc.hasAllele(genotypeAllele),
                                    "Allele " + genotypeAllele + " was present in genotype " + genotype +
                                            " but not in the VariantContext itself");
                        }
                    }
                }
            }
        }
    }

    /**
     * Helper method for testMaxAlternateAlleles
     *
     * @param vc VariantContext to check
     * @return number of alt alleles in vc, excluding NON_REF (if present)
     */
    private int getNumAltAllelesExcludingNonRef( final VariantContext vc ) {
        final List<Allele> altAlleles = vc.getAlternateAlleles();
        int numAltAllelesExcludingNonRef = 0;

        for ( final Allele altAllele : altAlleles ) {
            if ( ! altAllele.equals(Allele.NON_REF_ALLELE) ) {
                ++numAltAllelesExcludingNonRef;
            }
        }

        return numAltAllelesExcludingNonRef;
    }

    @Override
    public List<String> getToolSpecificArguments() {
        return Collections.emptyList();
    }
}
