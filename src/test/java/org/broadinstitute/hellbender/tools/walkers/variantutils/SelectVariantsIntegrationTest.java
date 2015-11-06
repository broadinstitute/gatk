package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

public class SelectVariantsIntegrationTest extends CommandLineProgramTest {

    private static String baseTestString(String args, String testFile) {
        return " --variant " + testFile
                    + " -o %s "
                    + args;
    }

    @Test
    public void testSampleSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + hg19MiniReference
                        + " --variant " + testFile
                        + " -sn NA11918 "
                        + " -sr " // suppress reference file name in output for test differencing
                        + " -o %s ",
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SimpleSelection.vcf")
        );

        spec.executeTest("testSampleSelection--" + testFile, this);
    }

    @Test
    public void testExpressionSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteringDepthInFormat.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + hg19MiniReference
                        + " --variant " + testFile
                        + " -select 'DP < 7' "
                        + " -sr " // suppress reference file name in output for test differencing
                        + " -o %s ",
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SimpleExpressionSelection.vcf")
        );

        spec.executeTest("testSimpleExpressionSelection--" + testFile, this);
    }

    @Test
    public void testRepeatedLineSelectionAndExludeFiltered() throws IOException {
        final String testFile = getToolTestDataDir() + "test.dup.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn A -sn B -sn C -ef ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_RepeatedLineSelection.vcf")
        );

        spec.executeTest("testRepeatedLineSelection--" + testFile, this);
    }

    @Test
    public void testComplexSelection()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.txt";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11894 -se 'NA069*' -sf " + samplesFile + " -select 'RMSMAPQ < 170.0'", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ComplexSelection.vcf")
        );

        spec.executeTest("testComplexSelection--" + testFile, this);
    }

    @Test
    public void testComplexSelectionWithNonExistingSamples()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.txt";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES  -select 'RMSMAPQ < 170.0' -sn Z -sf " // non existent samples on command line
                        + samplesFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ComplexSelectionWithNonExistingSamples.vcf")
        );
        spec.executeTest("testComplexSelectionWithNonExistingSamples--" + testFile, this);
    }

    @Test
    public void testNonExistingFieldSelection()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -env -select 'foo!=0 || RMSMAPQ < 170.0' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_NonExistingSelection.vcf")
        );

        spec.executeTest("testNonExistingSelection--" + testFile, this);
    }

    /**
     * Test excluding samples from file and sample name
     */
    @Test
    public void testSampleExclusionFromFileAndSeparateSample()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.txt";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xl_sn NA11894 -xl_sf " + samplesFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SampleExclusionFromFileAndSeparateSample.vcf")
        );

        spec.executeTest("testSampleExclusionFromFileAndSeparateSample--" + testFile, this);
    }

    /**
     * Test excluding samples from file
     */
    @Test
    public void testSampleExclusionJustFromFile()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.txt";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xl_sf " + samplesFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SampleExclusionJustFromFile.vcf")
        );

        spec.executeTest("testSampleExclusionJustFromFile--" + testFile, this);
    }

    /**
     * Test excluding samples from expression
     */
    @Test
    public void testSampleExclusionJustFromExpression()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xl_se 'NA069*' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SampleExclusionJustFromExpression.vcf")
        );

        spec.executeTest("testSampleExclusionJustFromExpression--" + testFile, this);
    }

    /**
     * Test excluding samples from negation expression
     */
    @Test
    public void testSampleExclusionJustFromNegationExpression()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -se 'NA[0-9]{4}[^1-9]' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SampleExclusionJustFromRegexExpression.vcf")
        );

        spec.executeTest("testSampleExclusionJustFromRegexExpression--" + testFile, this);
    }

    /**
     * Test including samples that are not in the VCF
     */

    @Test
    public void testSampleInclusionWithNonexistingSamples()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.txt";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn A -sn Z -sn Q -sf " + samplesFile, testFile),
                1,
                UserException.BadInput.class
        );

        spec.executeTest("testSampleInclusionWithNonexistingSamples--" + testFile, this);
    }

    @Test
    public void testDiscordance() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String discordanceFile = getToolTestDataDir() + "vcfexample2DiscordanceConcordance.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11992 " // not present in discordance track
                                + " -disc " + discordanceFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_Discordance.vcf")
        );

        spec.executeTest("testDiscordance--" + testFile, this);
    }

    @Test
    public void testConcordance()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String concordanceFile = getToolTestDataDir() + "vcfexample2DiscordanceConcordance.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11894 -conc " + concordanceFile  + " --lenient ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_Concordance.vcf")
        );

        spec.executeTest("testConcordance--" + testFile, this);
    }

    /**
     * Test including variant types.
     */
    @Test
    public void testVariantTypeSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -restrictAllelesTo MULTIALLELIC -selectType MIXED ",testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_VariantTypeSelection.vcf")
        );

        spec.executeTest("testVariantTypeSelection--" + testFile, this);
    }

    /**
     * Test excluding indels that are larger than the specified size
     */
    @Test
    public void testMaxIndelLengthSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -selectType INDEL --maxIndelSize 2 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MaxIndelLengthSelection.vcf")
        );

        spec.executeTest("testMaxIndelLengthSelection--" + testFile, this);
    }

    /**
     * Test excluding indels that are smaller than the specified size
     */
    @Test
    public void testMinIndelLengthSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -selectType INDEL --minIndelSize 2 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MinIndelLengthSelection.vcf")
        );

        spec.executeTest("testMinIndelLengthSelection--" + testFile, this);
    }

    @Test
    public void testRemoveMLE() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.withMLE.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA12892 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_RemoveMLE.vcf")
        );

        spec.executeTest("testRemoveMLE--" + testFile, this);
    }

    @Test
    public void testKeepOriginalAC() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.loseAlleleInSelection.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --keepOriginalAC -sn NA12892 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepOriginalAC.vcf")
        );

        spec.executeTest("testKeepOriginalAC--" + testFile, this);
    }

    @Test
    public void testKeepOriginalACAndENV() throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample.loseAlleleInSelection.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --keepOriginalAC -sn NA12892 -env -trimAlternates ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepOriginalACAndENV.vcf")
        );

        spec.executeTest("testKeepOriginalACAndENV--" + testFile, this);
    }

    @Test
    public void testKeepOriginalDP() throws IOException {
        final String testFile = getToolTestDataDir() + "CEUtrioTest.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --keepOriginalDP -sn NA12892 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepOriginalDP.vcf")
        );

        spec.executeTest("testKeepOriginalDP--" + testFile, this);
    }

    @Test
    public void testMultipleRecordsAtOnePosition() throws IOException {
        final String testFile = getToolTestDataDir() + "selectVariants.onePosition.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -select 'KG_FREQ < 0.5' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MultipleRecordsAtOnePosition.vcf")
        );

        spec.executeTest("testMultipleRecordsAtOnePosition--" + testFile, this);
    }

    @Test
    public void testNoGTs() throws IOException {
        final String testFile = getToolTestDataDir() + "vcf4.1.example.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec (
                " --variant " + testFile + " -o %s ",
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_NoGTs.vcf")
        );

        spec.executeTest("testNoGTs--" + testFile, this);
    }

    @Test
    public void testSelectFromMultiAllelic() throws IOException {
        final String testFile = getToolTestDataDir() + "multi-allelic.bi-allelicInGIH.vcf";
        final String samplesFile = getToolTestDataDir() + "GIH.samples.list";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sf " + samplesFile + " --excludeNonVariants -trimAlternates", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MultiAllelicExcludeNonVar.vcf")
        );
        spec.executeTest("test select from multi allelic with excludeNonVariants --" + testFile, this);
    }

    @Test
    public void testMultiAllelicAnnotationOrdering() throws IOException {
        final String testFile = getToolTestDataDir() + "multi-allelic-ordering.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn SAMPLE-CC -sn SAMPLE-CT -sn SAMPLE-CA --excludeNonVariants", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MultiAllelicAnnotationOrdering.vcf")
        );
        spec.executeTest("test multi allelic annotation ordering --" + testFile, this);
    }

    @Test
    public void testFileWithoutInfoLineInHeader() throws IOException {
        testFileWithoutInfoLineInHeader("testSelectVariants_FileWithoutInfoLineInHeader", IllegalStateException.class);
    }

    @Test
    public void testFileWithoutInfoLineInHeaderWithOverride() throws IOException {
        testFileWithoutInfoLineInHeader("testSelectVariants_FileWithoutInfoLineInHeaderWithOverride", null);
    }

    private void testFileWithoutInfoLineInHeader(final String name, final Class<? extends Exception> expectedException) throws IOException {
        final String testFile = getToolTestDataDir() + "missingHeaderLine.vcf";
        final String outFile = getToolTestDataDir() + "expected/" + name + ".vcf";

        final String cmd = baseTestString(" -sn NA12892 " + (expectedException == null ? " --lenient" : ""), testFile);

        IntegrationTestSpec spec =
                expectedException != null
                        ? new IntegrationTestSpec(cmd, 1, expectedException)
                        : new IntegrationTestSpec(cmd, Collections.singletonList(outFile));

        spec.executeTest(name, this);
    }

    @Test
    public void testInvalidJexl() throws IOException {
        final String testFile = getToolTestDataDir() + "ac0.vcf";

        // NOTE: JexlEngine singleton construction in VariantContextUtils sets silent to false.
        // However VariantFiltration.initialize() sets setSilent(true) on the shared instance.
        // Just in case this test runs after a VariantFiltration in the same VM, always set silent back to false.
        htsjdk.variant.variantcontext.VariantContextUtils.engine.get().setSilent(false);

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -select 'vc.getGenotype(\"FAKE_SAMPLE\").isHomRef()' ", testFile),
                1,
                UserException.class);
        spec.executeTest("InvalidJexl", this);
    }

    @Test
    public void testAlleleTrimming() throws IOException {
        final String testFile = getToolTestDataDir() + "forHardLeftAlignVariantsTest.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA12878 -env -trimAlternates ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_AlleleTrimming.vcf"));
        spec.executeTest("testAlleleTrimming", this);
    }

    @DataProvider(name="unusedAlleleTrimmingProvider")
    public Object[][] unusedAlleleTrimmingProvider() {
        final String expectedPath = getToolTestDataDir() + "expected/";
        return new Object[][] {
                {
                        getToolTestDataDir() + "forHardLeftAlignVariantsTest.vcf",
                        "-trimAlternates",
                        expectedPath + "testSelectVariants_UnusedAlleleHardLeftTrim.vcf"
                },
                {
                        getToolTestDataDir() + "forHardLeftAlignVariantsTest.vcf",
                        null,
                        expectedPath + "testSelectVariants_UnusedAlleleHardLeft.vcf"
                },
                {
                        getToolTestDataDir() + "multi-allelic-ordering.vcf",
                        "-sn SAMPLE-CC -sn SAMPLE-CT",
                        expectedPath + "testSelectVariants_UnusedAlleleCCCT.vcf"
                },
                {
                        getToolTestDataDir() + "multi-allelic-ordering.vcf",
                        "-sn SAMPLE-CC -sn SAMPLE-CT -env",
                        expectedPath + "testSelectVariants_UnusedAlleleCCCTEnv.vcf"
                },
                {
                        getToolTestDataDir() + "multi-allelic-ordering.vcf",
                        "-sn SAMPLE-CC -sn SAMPLE-CT -trimAlternates",
                        expectedPath + "testSelectVariants_UnusedAlleleCCCTTrim.vcf"
                },
                {
                        getToolTestDataDir() + "multi-allelic-ordering.vcf",
                        "-sn SAMPLE-CC -sn SAMPLE-CT -env -trimAlternates",
                        expectedPath + "testSelectVariants_UnusedAlleleCCCTTrimAltEnv.vcf"
                }
        };
    }

    @Test(dataProvider="unusedAlleleTrimmingProvider")
    public void testUnusedAlleleTrimming(final String vcf, final String extraArgs, final String expectedOutput) throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(extraArgs == null ? "" : extraArgs, vcf),
                Collections.singletonList(expectedOutput)
        );

        spec.executeTest(
                String.format("testUnusedAlleleTrimming: (%s,%s)", new File(vcf).getName(), extraArgs == null ? "(none)" : extraArgs),
                this);
    }

    /**
     *  Test with an empty VCF file
     */
    @Test
    public void testEmptyVcfException() throws IOException {
        final String testFile = getToolTestDataDir() + "reallyEmpty.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("", testFile),
                1,
                UserException.CouldNotReadInputFile.class
        );

        spec.executeTest("testEmptyVcfException--" + testFile, this);
    }

    /**
     * Test with a VCF file that is not a file
     */
    @Test
    public void testNotFileVcfException() throws IOException {
        final String testFile = getToolTestDataDir();

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("", testFile),
                1,
                UserException.CouldNotReadInputFile.class
        );

        spec.executeTest("testNotFileVcfException--" + testFile, this);
    }

    /**
     * Test with a VCF file that does not exist
     */
    @Test
    public void testMissingVcfException() throws IOException {
        final String testFile = "test.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("", testFile),
                1,
                UserException.CouldNotReadInputFile.class
        );

        spec.executeTest("testMissingVcfException--" + testFile, this);
    }

    /**
     * Test inverting the variant selection criteria by the -invertSelect argument
     */
    @Test
    public void testInvertSelection()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.txt";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11894 -sf " + samplesFile +
                                    " -select 'RMSMAPQ < 170.0' -invertSelect ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_InvertSelection.vcf")
        );

        spec.executeTest("testInvertSelection--" + testFile, this);
    }

    /**
     * Test inverting the variant selection criteria by inverting the JEXL expression logic following -select
     */
    @Test
    public void testInvertJexlSelection()  throws IOException {
        final String testFile = getToolTestDataDir() + "vcfexample2.vcf";
        final String samplesFile = getToolTestDataDir() + "samples.txt";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -sn NA11894 -sf " + samplesFile +
                        " -select 'RMSMAPQ > 170.0' ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_InvertJexlSelection.vcf")
        );

        spec.executeTest("testInvertJexlSelection--" + testFile, this);
    }

    /**
     * Test selecting variants with IDs
     */
    @Test
    public void testKeepSelectionID() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";
        final String idFile = getToolTestDataDir() + "complexExample1.vcf.id";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -IDs " + idFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_KeepSelectionID.vcf")
        );

        spec.executeTest("testKeepSelectionID--" + testFile, this);
    }

    /**
     * Test excluding variants with IDs
     */
    @Test
    public void testExcludeSelectionID() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";
        final String idFile = getToolTestDataDir() + "complexExample1.vcf.id";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xlIDs " + idFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ExcludeSelectionID.vcf")
        );

        spec.executeTest("testExcludeSelectionID--" + testFile, this);
    }

    /**
     * Test excluding variant types
     */
    @Test
    public void testExcludeSelectionType() throws IOException {
        final String testFile = getToolTestDataDir() + "complexExample1.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -xlSelectType SNP ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_ExcludeSelectionType.vcf")
        );

        spec.executeTest("testExcludeSelectionType--" + testFile, this);
    }

    @Test
    public void testMendelianViolationSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "CEUtrioTest.vcf";
        final String pedFile = getToolTestDataDir() + "CEUtrio.ped";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -ped " + pedFile + " -mv -mvq 0 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MendelianViolationSelection.vcf")
        );

        spec.executeTest("testMendelianViolationSelection--" + testFile, this);
    }

    @Test
    public void testInvertMendelianViolationSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "CEUtrioTest.vcf";
        final String pedFile = getToolTestDataDir() + "CEUtrio.ped";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" -mv -mvq 0 -invMv -ped " + pedFile, testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_InvertMendelianViolationSelection.vcf")
        );

        spec.executeTest("testInvertMendelianViolationSelection--" + testFile, this);
    }

    @Test
    public void testMaxFilteredGenotypesSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --maxFilteredGenotypes 1 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MaxFilteredGenotypesSelection.vcf")
        );

        spec.executeTest("testMaxFilteredGenotypesSelection--" + testFile, this);
    }

    @Test
    public void testMinFilteredGenotypesSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --minFilteredGenotypes 2 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MinFilteredGenotypesSelection.vcf")
        );

        spec.executeTest("testMinFilteredGenotypesSelection--" + testFile, this);
    }

    @Test
    public void testMaxFractionFilteredGenotypesSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --maxFractionFilteredGenotypes 0.4 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MaxFractionFilteredGenotypesSelection.vcf")
        );

        spec.executeTest("testMaxFractionFilteredGenotypesSelection--" + testFile, this);
    }

    @Test
    public void testMinFractionFilteredGenotypesSelection() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final  IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --minFractionFilteredGenotypes 0.6 ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_MinFractionFilteredGenotypesSelection.vcf")
        );

        spec.executeTest("testMinFractionFilteredGenotypesSelection--" + testFile, this);
    }

    @Test
    public void testSetFilteredGtoNocall() throws IOException {
        final String testFile = getToolTestDataDir() + "filteredSamples.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(" --setFilteredGtToNocall ", testFile),
                Collections.singletonList(getToolTestDataDir() + "expected/" + "testSelectVariants_SetFilteredGtoNocall.vcf")
        );

        spec.executeTest("testSetFilteredGtoNocall--" + testFile, this);
    }
}
