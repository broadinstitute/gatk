package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SamFiles;
import htsjdk.tribble.Tribble;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class HaplotypeCallerIntegrationTest extends CommandLineProgramTest {

    public static final String TEST_FILES_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/haplotypecaller/";

    @DataProvider(name="HaplotypeCallerTestInputs")
    public Object[][] getHaplotypCallerTestInputs() {
        return new Object[][] {
                {NA12878_20_21_WGS_bam, b37_reference_20_21},
                {NA12878_20_21_WGS_cram, b37_reference_20_21}
        };
    }
    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testVCFModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testVCFMode.gatk4.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "-addOutputVCFCommandLine", "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }


    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     *
     * Test currently throws an exception due to lack of support for allele-specific annotations in VCF mode
     */
    @Test(dataProvider="HaplotypeCallerTestInputs", expectedExceptions = UserException.class)
    public void testVCFModeIsConsistentWithPastResults_AlleleSpecificAnnotations(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        //NOTE: AlleleSpecific support in the VCF mode is implemented but bogus for now.
        // This test should not be treated as a strict check of correctness.
        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk4.alleleSpecific.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "-addOutputVCFCommandLine", "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    /*
     * Test that in VCF mode we're >= 99% concordant with GATK3.8 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testVCFModeIsConcordantWithGATK3_8Results(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3Results", ".vcf");
        //Created by running:
        // java -jar gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // --out expected.testVCFMode.gatk3.8-4-g7b0250253.vcf -G StandardHC -G Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk3.8-4-g7b0250253.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in VCF mode is < 99% (" +  concordance + ")");
    }

    /*
     * Test that in VCF mode we're >= 99% concordant with GATK3.8 results
     *
     * Test currently throws an exception due to lack of support for allele-specific annotations in VCF mode
     */
    @Test(dataProvider="HaplotypeCallerTestInputs", expectedExceptions = UserException.class)
    public void testVCFModeIsConcordantWithGATK3_8ResultsAlleleSpecificAnnotations(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3.8ResultsAlleleSpecificAnnotations", ".vcf");

        //Created by running
        //java -jar gatk.3.8-4-g7b0250253f.jar  -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // --out expected.testVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.vcf -G StandardHC -G Standard -G AS_Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "-pairHMM", "AVX_LOGLESS_CACHING",
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in AS VCF mode is < 99% (" +  concordance + ")");
    }

    /*
     * Test that in GVCF mode we're consistent with past GATK4 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testGVCFModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults", ".g.vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testGVCFMode.gatk4.g.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "-addOutputVCFCommandLine", "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
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

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "-addOutputVCFCommandLine", "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    /*
     * Test that in GVCF mode we're >= 99% concordant with GATK3 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testGVCFModeIsConcordantWithGATK3_8Results(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3Results", ".g.vcf");
        //Created by running:
        //java -jar  gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // -ERC GVCF --out expected.testGVCFMode.3.8-4-g7b0250253f.g.vcf -G StandardHC -G Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testGVCFMode.3.8-4-g7b0250253f.g.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING",
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in GVCF mode is < 99% (" +  concordance + ")");
    }

    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testGVCFModeIsConcordantWithGATK3_8AlelleSpecificResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();
        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3_8AlelleSpecificResults", ".g.vcf");

        //Created by running:
        // java -jar gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // -ERC GVCF --out expected.testGVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.g.vcf -G StandardHC -G Standard -G AS_Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testGVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.g.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING",
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in AS GVCF mode is < 99% (" +  concordance + ")");
    }

    @Test
    public void testBamoutProducesReasonablySizedOutput() {
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
        final File bamOutput = createTempFile("testBamoutProducesReasonablySizedOutput", ".bam");

        ArgumentsBuilder argBuilder = new ArgumentsBuilder();

        argBuilder.addInput(new File(NA12878_20_21_WGS_bam));
        argBuilder.addReference(new File(b37_reference_20_21));
        argBuilder.addOutput(new File(vcfOutput.getAbsolutePath()));
        argBuilder.addArgument("L", testInterval);
        argBuilder.addArgument("bamout", bamOutput.getAbsolutePath());
        argBuilder.addArgument("pairHMM", "AVX_LOGLESS_CACHING");

        runCommandLine(argBuilder.getArgsArray());

        try ( final ReadsDataSource bamOutReadsSource = new ReadsDataSource(bamOutput.toPath()) ) {
            int actualBamoutNumReads = 0;
            for ( final GATKRead read : bamOutReadsSource ) {
                ++actualBamoutNumReads;
            }

            final int readCountDifference = Math.abs(actualBamoutNumReads - gatk3BamoutNumReads);
            Assert.assertTrue(((double)readCountDifference / gatk3BamoutNumReads) < 0.10,
                    "-bamout produced a bam with over 10% fewer/more reads than expected");
        }
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
            final boolean createVCFOutMD5) {
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
        argBuilder.addArgument("bamout", bamOutput.getAbsolutePath());
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
    }

    @Test
    public void testHaplotypeCallerRemoveAltAlleleBasedOnHaptypeScores() throws IOException{
        final File testBAM = new File(TEST_FILES_DIR + "pretendTobeTetraPloidTetraAllelicSite.bam");
        final File output = createTempFile("testHaplotypeCallerRemoveAltAlleleBasedOnHaptypeScoresResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testHaplotypeCallerRemoveAltAlleleBasedOnHaptypeScores.gatk4.vcf");

        final String[] args = {
                "-I", testBAM.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20:11363580-11363600",
                "-O", output.getAbsolutePath(),
                "-ploidy", "4",
                "-maxGT", "15",
                "-addOutputVCFCommandLine", "false"
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    // test that ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT removes reads that consume zero reference bases
    // e.g. read name HAVCYADXX150109:1:2102:20528:2129 with cigar 23S53I
    @Test
    public void testReadsThatConsumeZeroReferenceReads() throws Exception {
        final String CONSUMES_ZERO_REFERENCE_BASES = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/na12878-chr20-consumes-zero-reference-bases.bam";
        final File outputVcf = createTempFile("output", ".vcf");
        final String[] args = {
                "-I", CONSUMES_ZERO_REFERENCE_BASES,
                "-R", b37_reference_20_21,
                "-O", outputVcf.getAbsolutePath()
        };
        runCommandLine(args);
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
                "--assemblyRegionOut", assemblyRegionOut.getAbsolutePath(),
                "--activityProfileOut", activityProfileOut.getAbsolutePath()
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(assemblyRegionOut, expectedAssemblyRegionOut);
        IntegrationTestSpec.assertEqualTextFiles(activityProfileOut, expectedActivityProfileOut);
    }

    /*
     * Calculate rough concordance between two vcfs, comparing only the positions, alleles, and the first genotype.
     */
    public static double calculateConcordance( final File actual, final File expected ) {
        final Set<String> actualVCFKeys = new HashSet<>();
        final Set<String> expectedVCFKeys = new HashSet<>();
        int concordant = 0;
        int discordant = 0;

        try ( final FeatureDataSource<VariantContext> actualSource = new FeatureDataSource<>(actual);
              final FeatureDataSource<VariantContext> expectedSource = new FeatureDataSource<>(expected) ) {

            for ( final VariantContext vc : actualSource ) {
                actualVCFKeys.add(keyForVariant(vc));
            }

            for ( final VariantContext vc : expectedSource ) {
                expectedVCFKeys.add(keyForVariant(vc));
            }

            for ( final String vcKey : actualVCFKeys ) {
                if ( ! expectedVCFKeys.contains(vcKey) ) {
                    ++discordant;
                }
                else {
                    ++concordant;
                }
            }

            for ( final String vcKey : expectedVCFKeys ) {
                if ( ! actualVCFKeys.contains(vcKey) ) {
                    ++discordant;
                }
            }
        }

        return (double)concordant / (double)(concordant + discordant);
    }

    private static String keyForVariant( final VariantContext variant ) {
        return String.format("%s:%d-%d %s %s", variant.getContig(), variant.getStart(), variant.getEnd(),
                variant.getAlleles(), variant.getGenotype(0).getGenotypeString(false));
    }
}
