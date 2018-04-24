package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SamFiles;
import htsjdk.tribble.Tribble;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class HaplotypeCallerIntegrationTest extends CommandLineProgramTest {

    public static final String TEST_FILES_DIR = toolsTestDir + "haplotypecaller/";

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
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
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
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
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
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
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
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
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
    public void testGenotypeGivenAllelesMode() throws IOException {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGenotypeGivenAllelesMode", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testGenotypeGivenAllelesMode.gatk4.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                "--genotyping-mode", "GENOTYPE_GIVEN_ALLELES",
                "--alleles", new File(TEST_FILES_DIR, "testGenotypeGivenAllelesMode_givenAlleles.vcf").getAbsolutePath()
        };

        runCommandLine(args);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
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
        argBuilder.addArgument(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_SHORT_NAME, bamOutput.getAbsolutePath());
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
                "--max-genotype-count", "15",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
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

    // Test fix for https://github.com/broadinstitute/gatk/issues/3466
    // This specifically tests the case of a read that ends up as empty
    // after a call to ReadClipper.hardClipSoftClippedBases()
    @Test
    public void testCompletelyClippedReadNearStartOfContig() throws Exception {
        final File testCaseFilesDir = new File(TEST_FILES_DIR, "issue3466_gatk_cigar_error");
        final File output = createTempFile("testCompletelyClippedReadNearStartOfContig", ".vcf");
        final File expected = new File(testCaseFilesDir, "expected_output_gatk3.vcf");

        final String[] args = {
                "-I", new File(testCaseFilesDir, "culprit.bam").getAbsolutePath(),
                "-R", new File(testCaseFilesDir, "GRCh37_MTonly.fa").getAbsolutePath(),
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);

        Assert.assertEquals(calculateConcordance(output, expected), 1.0);
    }

    // Test fix for https://github.com/broadinstitute/gatk/issues/3845
    // This specifically tests the case of a read at the start of a contig (position == 1)
    // that becomes completely clipped after a call to ReadClipper.revertSoftClippedBases()
    @Test
    public void testCompletelyClippedReadNearStartOfContig_revertSoftClipped() throws Exception {
        final File testCaseFilesDir = new File(TEST_FILES_DIR, "issue3845_revertSoftClip_bug");
        final File output = createTempFile("testCompletelyClippedReadNearStartOfContig_revertSoftClipped", ".vcf");
        final File expected = new File(testCaseFilesDir, "expected_testCompletelyClippedReadNearStartOfContig_revertSoftClipped_gatk4.vcf");

        final String[] args = {
                "-I", new File(testCaseFilesDir, "issue3845_bug.bam").getAbsolutePath(),
                "-R", new File(publicTestDir, "Homo_sapiens_assembly38_chrM_only.fasta").getAbsolutePath(),
                "-L", "chrM",
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);

        Assert.assertEquals(calculateConcordance(output, expected), 1.0);
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

    /*
    * Test that the min_base_quality_score parameter works
    */
    @Test
    public void testMinBaseQualityScore() throws Exception {
        Utils.resetRandomGenerator();

        final File outputAtLowThreshold = createTempFile("output", ".vcf");
        final File outputAtHighThreshold = createTempFile("output", ".vcf");

        final String[] lowThresholdArgs = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", outputAtLowThreshold.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME, "20"
        };

        runCommandLine(lowThresholdArgs);

        final String[] highThresholdArgs = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", outputAtHighThreshold.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME, "30"
        };

        runCommandLine(highThresholdArgs);

        try (final FeatureDataSource<VariantContext> lowThresholdSource = new FeatureDataSource<VariantContext>(outputAtLowThreshold);
             final FeatureDataSource<VariantContext> highThresholdSource = new FeatureDataSource<VariantContext>(outputAtHighThreshold)) {
            final List<VariantContext> variantsWithLowThreshold =
                    StreamSupport.stream(lowThresholdSource.spliterator(), false).collect(Collectors.toList());

            final List<VariantContext> variantsWithHighThreshold =
                    StreamSupport.stream(highThresholdSource.spliterator(), false).collect(Collectors.toList());

            final Set<Integer> lowStarts = variantsWithLowThreshold.stream().map(VariantContext::getStart).collect(Collectors.toSet());
            final Set<Integer> highStarts = variantsWithHighThreshold.stream().map(VariantContext::getStart).collect(Collectors.toSet());
            lowStarts.removeAll(highStarts);
            final List<Integer> diff = lowStarts.stream().sorted().collect(Collectors.toList());
            Assert.assertEquals(diff, Arrays.asList(10002458, 10008758, 10009842));
        }

    }

    @DataProvider
    public Object[][] getContaminationCorrectionTestData() {
        // This bam is a snippet of GATKBaseTest.NA12878_20_21_WGS_bam artificially contaminated
        // at 15% with reads from another sample. We use such a high contamination fraction to ensure
        // that calls will be noticeably different vs. running without -contamination at all.
        final String contaminatedBam15Percent = largeFileTestDir + "contaminated_bams/CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.bam";

        final SimpleInterval traversalInterval = new SimpleInterval("20", 10100000, 10150000);

        // File specifying the per-sample contamination for contaminatedBam15Percent
        final String contaminationFile = TEST_FILES_DIR + "contamination_file_for_CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT";

        // Expected calls from GATK 3.x on the contaminated bam running with -contamination
        // Created in GATK 3.8-1-1-gdde23f56a6 using the command:
        // java -jar target/GenomeAnalysisTK.jar -T HaplotypeCaller -I ../hellbender/src/test/resources/large/contaminated_bams/CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.bam -R ../hellbender/src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10100000-10150000 -contamination 0.15 -o ../hellbender/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.calls.20.10100000-10150000.gatk3.8-1-1-gdde23f56a6.vcf
        final String expectedGATK3ContaminationCorrectedCallsVCF = TEST_FILES_DIR + "expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.calls.20.10100000-10150000.gatk3.8-1-1-gdde23f56a6.vcf";

        // Expected calls from GATK 3.x on the contaminated bam running with -contamination
        // Created in GATK 3.8-1-1-gdde23f56a6 using the command:
        // java -jar target/GenomeAnalysisTK.jar -T HaplotypeCaller -I ../hellbender/src/test/resources/large/contaminated_bams/CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.bam -R ../hellbender/src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10100000-10150000 -contamination 0.15 -o ../hellbender/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.calls.20.10100000-10150000.gatk3.8-1-1-gdde23f56a6.g.vcf -ERC GVCF
        final String expectedGATK3ContaminationCorrectedCallsGVCF = TEST_FILES_DIR + "expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.calls.20.10100000-10150000.gatk3.8-1-1-gdde23f56a6.g.vcf";

        // Expected calls from GATK 4 on the uncontaminated bam (VCF mode)
        final String expectedGATK4UncontaminatedCallsVCF = TEST_FILES_DIR + "expected.CEUTrio.HiSeq.WGS.b37.NA12878.calls.20.10100000-10150000.vcf";

        // Expected calls from GATK 4 on the uncontaminated bam (GVCF mode)
        final String expectedGATK4UncontaminatedCallsGVCF = TEST_FILES_DIR + "expected.CEUTrio.HiSeq.WGS.b37.NA12878.calls.20.10100000-10150000.g.vcf";

        // bam,
        // known contamination fraction,
        // per-sample contamination file,
        // interval,
        // reference,
        // GVCF mode,
        // GATK 4 calls on the uncontaminated bam
        // GATK 3 calls on the contaminated bam with contamination correction
        return new Object[][] {
                { contaminatedBam15Percent,
                  0.15,
                  contaminationFile,
                  traversalInterval,
                  b37_reference_20_21,
                  false, // VCF mode
                  expectedGATK4UncontaminatedCallsVCF,
                  expectedGATK3ContaminationCorrectedCallsVCF
                },
                { contaminatedBam15Percent,
                  0.15,
                  contaminationFile,
                  traversalInterval,
                  b37_reference_20_21,
                  true, // GVCF mode
                  expectedGATK4UncontaminatedCallsGVCF,
                  expectedGATK3ContaminationCorrectedCallsGVCF
                }
        };
    }

    @Test(dataProvider = "getContaminationCorrectionTestData")
    public void testContaminationCorrection( final String contaminatedBam,
                                   final double contaminationFraction,
                                   final String contaminationFile,
                                   final SimpleInterval interval,
                                   final String reference,
                                   final boolean gvcfMode,
                                   final String gatk4UncontaminatedCallsVCF,
                                   final String gatk3ContaminationCorrectedCallsVCF ) throws Exception {
        final File uncorrectedOutput = createTempFile("testContaminationCorrectionUncorrectedOutput", gvcfMode ? ".g.vcf" : ".vcf");
        final File correctedOutput = createTempFile("testContaminationCorrectionCorrectedOutput", gvcfMode ? ".g.vcf" : ".vcf");
        final File correctedOutputUsingContaminationFile = createTempFile("testContaminationCorrectionCorrectedOutputUsingContaminationFile", gvcfMode ? ".g.vcf" : ".vcf");

        // Generate raw uncorrected calls on the contaminated bam, for comparison purposes
        final String[] noContaminationCorrectionArgs = {
                "-I", contaminatedBam,
                "-R", reference,
                "-L", interval.toString(),
                "-O", uncorrectedOutput.getAbsolutePath(),
                "-ERC", (gvcfMode ? "GVCF" : "NONE")
        };
        Utils.resetRandomGenerator();
        runCommandLine(noContaminationCorrectionArgs);

        // Run with contamination correction using -contamination directly
        final String[] contaminationCorrectionArgs = {
                "-I", contaminatedBam,
                "-R", reference,
                "-L", interval.toString(),
                "-O", correctedOutput.getAbsolutePath(),
                "-contamination", Double.toString(contaminationFraction),
                "-ERC", (gvcfMode ? "GVCF" : "NONE")
        };
        Utils.resetRandomGenerator();
        runCommandLine(contaminationCorrectionArgs);

        // Run with contamination correction using a contamination file
        final String[] contaminationCorrectionFromFileArgs = {
                "-I", contaminatedBam,
                "-R", reference,
                "-L", interval.toString(),
                "-O", correctedOutputUsingContaminationFile.getAbsolutePath(),
                "-contamination-file", contaminationFile,
                "-ERC", (gvcfMode ? "GVCF" : "NONE")
        };
        Utils.resetRandomGenerator();
        runCommandLine(contaminationCorrectionFromFileArgs);

        // Calculate concordance vs. the calls on the uncontaminated bam. Count only actual calls
        // (ignoring GVCF reference blocks, if present) for the purposes of the concordance calculation.
        final double uncorrectedCallsVSUncontaminatedCallsConcordance = calculateConcordance(uncorrectedOutput, new File(gatk4UncontaminatedCallsVCF), gvcfMode);
        final double correctedCallsVSUncontaminatedCallsConcordance = calculateConcordance(correctedOutput, new File(gatk4UncontaminatedCallsVCF), gvcfMode);
        final double correctedCallsFromFileVSUncontaminatedCallsConcordance = calculateConcordance(correctedOutputUsingContaminationFile, new File(gatk4UncontaminatedCallsVCF), gvcfMode);

        // Calculate concordance vs. the contamination-corrected calls from GATK3. Here we count
        // all records in the output, including reference blocks.
        final double correctedCallsVSGATK3CorrectedCallsConcordance = calculateConcordance(correctedOutput, new File(gatk3ContaminationCorrectedCallsVCF));
        final double correctedCallsFromFileVSGATK3CorrectedCallsConcordance = calculateConcordance(correctedOutputUsingContaminationFile, new File(gatk3ContaminationCorrectedCallsVCF));

        // Sanity checks: the concordance when running with -contamination should be the same as the
        // concordance when running with -contamination-file
        assertEqualsDoubleSmart(correctedCallsVSUncontaminatedCallsConcordance, correctedCallsFromFileVSUncontaminatedCallsConcordance, 0.001,
                "concordance running with -contamination and -contamination-file should be identical, but wasn't");
        assertEqualsDoubleSmart(correctedCallsVSGATK3CorrectedCallsConcordance, correctedCallsFromFileVSGATK3CorrectedCallsConcordance, 0.001,
                "concordance running with -contamination and -contamination-file should be identical, but wasn't");

        // With -contamination correction on, concordance vs. the calls on the uncontaminated bam
        // should improve by at least 1% (actual improvement tends to be much larger, but we just
        // want to make sure that the calls don't actually get worse here)
        Assert.assertTrue(correctedCallsVSUncontaminatedCallsConcordance > uncorrectedCallsVSUncontaminatedCallsConcordance + 0.01,
                "running with -contamination should have improved concordance vs. calls on the uncontaminated bam by at least 1%, but it didn't");

        // With -contamination-file correction on, concordance vs. the calls on the uncontaminated bam
        // should improve by at least 1% (actual improvement tends to be much larger, but we just
        // want to make sure that the calls don't actually get worse here)
        Assert.assertTrue(correctedCallsFromFileVSUncontaminatedCallsConcordance > uncorrectedCallsVSUncontaminatedCallsConcordance + 0.01,
                "running with -contamination-file should have improved concordance vs. calls on the uncontaminated bam by at least 1%, but it didn't");

        // -contamination corrected output in GATK4 should have >= 99% concordance
        // vs. -contamination corrected output in GATK3
        Assert.assertTrue(correctedCallsVSGATK3CorrectedCallsConcordance >= 0.99,
                "output with -contamination correction should have >= 99% concordance with contamination-corrected calls from GATK3, but it didn't");

        // -contamination-file corrected output in GATK4 should have >= 99% concordance
        // vs. -contamination corrected output in GATK3
        Assert.assertTrue(correctedCallsFromFileVSGATK3CorrectedCallsConcordance >= 0.99,
                "output with -contamination-file correction should have >= 99% concordance with contamination-corrected calls from GATK3, but it didn't");
    }

    // With 100% contamination in VCF mode, make sure that we get no calls
    @Test
    public void test100PercentContaminationNoCallsInVCFMode() throws Exception {
        final File output = createTempFile("test100PercentContaminationNoCallsInVCFMode", ".vcf");

        final String[] contaminationArgs = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", output.getAbsolutePath(),
                "-contamination", "1.0"
        };
        runCommandLine(contaminationArgs);

        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        // Check that we get no calls, but a non-empty VCF header.
        Assert.assertTrue(result.getLeft().getMetaDataInInputOrder().size() > 0, "There should be a non-empty header present");
        Assert.assertEquals(result.getRight().size(), 0, "There should be no calls with 100% contamination in VCF mode");
    }

    // With 100% contamination in GVCF mode, make sure that we get no calls (only reference blocks)
    @Test
    public void test100PercentContaminationNoCallsInGVCFMode() throws Exception {
        final File output = createTempFile("test100PercentContaminationNoCallsInGVCFMode", ".g.vcf");

        final String[] contaminationArgs = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", output.getAbsolutePath(),
                "-contamination", "1.0",
                "-ERC", "GVCF"
        };
        runCommandLine(contaminationArgs);

        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        // Check that we get a non-empty VCF header.
        Assert.assertTrue(result.getLeft().getMetaDataInInputOrder().size() > 0, "There should be a non-empty header present");

        // Check that we get only GVCF reference blocks, and no actual calls
        final List<VariantContext> vcfRecords = result.getRight();
        Assert.assertTrue(! vcfRecords.isEmpty(), "VCF should be non-empty in GVCF mode");
        for ( final VariantContext vc : vcfRecords ) {
            Assert.assertTrue(isGVCFReferenceBlock(vc), "Expected only GVCF reference blocks (no actual calls)");
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

    /*
     * Calculate rough concordance between two vcfs, comparing only the positions, alleles, and the first genotype.
     */
    public static double calculateConcordance( final File actual, final File expected ) {
        return calculateConcordance(actual, expected, false);
    }

    /*
     * Calculate rough concordance between two vcfs, comparing only the positions, alleles, and the first genotype.
     * Can optionally ignore GVCF blocks in the concordance calculation.
     */
    public static double calculateConcordance( final File actual, final File expected, final boolean ignoreGVCFBlocks ) {
        final Set<String> actualVCFKeys = new HashSet<>();
        final Set<String> expectedVCFKeys = new HashSet<>();
        int concordant = 0;
        int discordant = 0;

        try ( final FeatureDataSource<VariantContext> actualSource = new FeatureDataSource<>(actual);
              final FeatureDataSource<VariantContext> expectedSource = new FeatureDataSource<>(expected) ) {

            for ( final VariantContext vc : actualSource ) {
                if ( ! ignoreGVCFBlocks || ! isGVCFReferenceBlock(vc) ) {
                    actualVCFKeys.add(keyForVariant(vc));
                }
            }

            for ( final VariantContext vc : expectedSource ) {
                if ( ! ignoreGVCFBlocks || ! isGVCFReferenceBlock(vc) ) {
                    expectedVCFKeys.add(keyForVariant(vc));
                }
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

    private static boolean isGVCFReferenceBlock( final VariantContext vc ) {
        return vc.hasAttribute(VCFConstants.END_KEY) &&
               vc.getAlternateAlleles().size() == 1 &&
               vc.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE);
    }
}
