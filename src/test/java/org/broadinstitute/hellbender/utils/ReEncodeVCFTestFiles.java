package org.broadinstitute.hellbender.utils;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

// Issues that surface when fully decoding without using lenient processing:
//
// src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/tetra-diploid.vcf:
// htsjdk.tribble.TribbleException$InvalidHeader: Your input file has a malformed header: Discordant field size
// detected for field MLEAC at 20:10316239.  Field had 1 values but the header says this should have 2 values based
// on header record INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the
// allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
//
//src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testMask4.vcf:
// htsjdk.tribble.TribbleException$InvalidHeader: Your input file has a malformed header: Discordant field size
// detected for field AS_FilterStatus at 1:10020416.  Field had 1 values but the header says this should have 2
// values based on header record INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description="Filter status for each
// allele, as assessed by ApplyVQSR. Note that the VCF filter field will reflect the most lenient/sensitive status
// across all alleles.">
//
//src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_split_with_AS_filters.vcf
// htsjdk.tribble.TribbleException$InvalidHeader: Your input file has a malformed header: Discordant field size
// detected for field AS_SB_TABLE at chrM:73.  Field had 3 values but the header says this should have 1 values based
// on header record INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read
// counts for strand bias tests">

public class ReEncodeVCFTestFiles extends CommandLineProgramTest {

    @Override
    public String getTestedToolName() { return "SelectVariants"; }

    @DataProvider
    public Object[][] vcfFilesToReEncode() {
        return new Object[][]{
                //NOTE: In addition to these files, there is one remote test file (testSelectVariants_SimpleSelection.vcf)
                // in the gcp staging area that will need to be updated.
                //
                // Test input files that need to be re-encoded to get outputs to match re-encoded outputs
                // Some of these files could not be re-encoded (or rather, fully decoded) without modifications, it., to the header
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/CalculateGenotypePosteriors/NA12878.Jan2013.haplotypeCaller.subset.indels.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/CalculateGenotypePosteriors/CEUtrioTest_chr1.vcf"},
//                {"src/test/resources/large/1000G.phase3.broad.withGenotypes.chr20.10100000.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/CalculateGenotypePosteriors/overlappingVariants.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/test_left_align_hg38.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/gvcfForReblocking.g.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/expected.aggressiveQualFiltering.g.vcf"},
//                // required addition of MLEAC/MLEAF info header lines in order to fully decode
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/nonRefAD.g.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/testOneSampleAsForGnomAD.expected.g.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/prod.chr20snippet.withRawMQ.g.vcf"},

                // Expected output files that will get re-encoded as part of the 4.3 upgrade.
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/vcfexample2.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/complexExample1.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/haploid-multisample.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/CEUtrioTest.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/vcfexample.forNoCallFiltering.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/multi-allelic-ordering.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/selectVariants.onePosition.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/vcf4.1.example.vcf"},
//                // note:this is v3.3 file and has to be re-encoded using SelectVariants
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/test.dup.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/261_S01_raw_variants_gvcf.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/multi-allelic.bi-allelicInGIH.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/tetra-diploid.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/tetraploid-multisample.vcf"},
//                {"src/test/resources/large/CalculateGenotypePosteriors/expectedCGP_testMissingPriors.vcf"},
//                {"src/test/resources/large/CalculateGenotypePosteriors/expectedCGP_testUsingDiscoveredAF.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testClusteredSnps.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testFilter1.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testFilter2.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testFilterWithSeparateNames.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testGenotypeFilters1.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testGenotypeFilters2.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testInvertFilter.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testInvertJexlFilter.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testMask1.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testMask2.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testMask3.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testMask4.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testMaskReversed.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration/expected/testVariantFiltration_testNoAction.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/CalculateGenotypePosteriors/expectedCGP_testDefaultsWithPanel.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/CalculateGenotypePosteriors/expectedCGP_testFamilyPriors_chr1.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/CalculateGenotypePosteriors/expectedCGP_testInputINDELs.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/CalculateGenotypePosteriors/expectedCGP_testNumRefWithPanel.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/CalculateGenotypePosteriors/overlappingVariants.expected.no_duplicates.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_left_align_hg38.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_left_align_hg38_maxIndelSize296.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_left_align_hg38_maxIndelSize342.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_left_align_hg38_notrim.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_left_align_hg38_notrim_split_multiallelics.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_left_align_hg38_split_multiallelics.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_left_align_hg38_split_multiallelics_keepOrigAC.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/LeftAlignAndTrimVariants/expected_split_with_AS_filters.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/prod.chr20snippet.withRawMQ.expected.g.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/testJustOneSample.expected.g.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/testNonRefADCorrection.expected.g.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_ComplexSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_ComplexSelectionWithNonExistingSamples.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_Concordance.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_Discordance.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_DropAnnotations.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_DropAnnotationsSelectFisherStrand.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_DropAnnotationsSelectGQ.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_DropAnnotationsSelectRD.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_DropAnnotationsSelectRMSMAPQ.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_ExcludeSelectionID.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_ExcludeSelectionType.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_Haploid.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_InvertJexlSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_InvertMendelianViolationSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_InvertSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_KeepOriginalDP.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_KeepSelectionID.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_MaxFilteredGenotypesSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_MendelianViolationSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_MinIndelLengthSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_MultiAllelicAnnotationOrdering.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_MultiAllelicExcludeNonVar.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_MultipleRecordsAtOnePosition.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_NoGTs.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_RemoveSingleSpanDelAlleleExNoVar.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_RepeatedLineSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_SampleExclusionFromFileAndSeparateSample.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_SampleExclusionJustFromExpression.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_SampleExclusionJustFromFile.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_SampleExclusionJustFromRegexExpression.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_SimpleDiploid.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_SimpleSelection.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_TetraDiploid.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_Tetraploid.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_UnusedAlleleCCCT.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_UnusedAlleleCCCTEnv.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_UnusedAlleleCCCTTrim.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_UnusedAlleleCCCTTrimAltEnv.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_maxNOCALLnumber1.vcf"},
//                {"src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/SelectVariants/expected/testSelectVariants_maxNOCALLnumber2.vcf"},
        };
    }

    // This is not really a test, just a convenient way to house code. Reencodes the VCF files in the data
    // provider using SelectVariants.
    // NOTE: this generates an index file whether one existed before or not
    @Test(dataProvider = "vcfFilesToReEncode"
            ,enabled = false
    ) // don't actually run this test during CI
    public void reencodeVCFFile(final String sourceVCFName) throws IOException {
        final GATKPath sourceVCFPath = new GATKPath(sourceVCFName);
        final String tempFileName = sourceVCFPath.toPath().toFile().getName();
        final File tempVCFFile = IOUtils.createTempFile(tempFileName, (new GATKPath(sourceVCFName).getExtension().get()));

        // Copy the original file to use as the input, and re-encoded by rewriting the original.
        FileUtils.copyFile(sourceVCFPath.toPath().toFile(), tempVCFFile, false);

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addVCF(tempVCFFile);
        argBuilder.addOutput(sourceVCFPath.toPath());
        argBuilder.add(StandardArgumentDefinitions.LENIENT_LONG_NAME, true);
        argBuilder.add("fully-decode", true);
        argBuilder.add(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME, true);
        argBuilder.add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(argBuilder.getArgsArray());
    }
}