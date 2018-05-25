package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.mafOutput.MafOutputRendererConstants;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription;

/**
 * An integration test for the {@link Funcotator} tool.
 * Created by jonn on 8/29/17.
 */
public class FuncotatorIntegrationTest extends CommandLineProgramTest {

    //TODO: SUPER IMPORTANT! Add checks on the output file to make sure it's equal to an expected output file!

    // Temp directory in which to place output files.
    private static final File tmpOutDir;

    // Whether to do debug output (i.e. leave output around).
    // This should always be false when checked in.
    private static final boolean doDebugTests             = false;
    private static final String  LARGE_DATASOURCES_FOLDER = "funcotator_dataSources_latest";

    private static final String PIK3CA_VCF_HG19 = toolsTestDir + "funcotator/0816201804HC0_R01C01.pik3ca.vcf";
    private static final String PIK3CA_VCF_HG38 = toolsTestDir + "funcotator/hg38_trio.pik3ca.vcf";
    private static final String PIK3CA_VCF_HG19_SNPS = toolsTestDir + "funcotator/PIK3CA_SNPS_3.vcf";
    private static final String PIK3CA_VCF_HG19_INDELS = toolsTestDir + "funcotator/PIK3CA_INDELS_3.vcf";
    private static final String MUC16_VCF_HG19 = toolsTestDir + "funcotator/MUC16_MNP.vcf";
    private static final String PIK3CA_VCF_HG19_ALTS = toolsTestDir + "funcotator/PIK3CA_3_miss_clinvar_alt_only.vcf";
    private static final String DS_PIK3CA_DIR = largeFileTestDir + "funcotator/small_ds_pik3ca/";
    private static final String DS_MUC16_DIR = largeFileTestDir + "funcotator/small_ds_muc16/";

    private static String hg38Chr3Ref;
    private static String b37Chr3Ref;
    private static String hg19Chr3Ref;
    private static String hg19Chr19Ref;
    static {
        if (!doDebugTests) {
            tmpOutDir = createTempDir("funcotatorTmpFolder");
        } else {
            tmpOutDir = new File("funcotatorTmpFolder" + File.separator);
            if (!tmpOutDir.mkdirs() && !tmpOutDir.exists()) {
                throw new GATKException("Error making output folder for test: " + FuncotatorIntegrationTest.class.getName());
            }
        }

        hg38Chr3Ref = FuncotatorReferenceTestUtils.retrieveHg38Chr3Ref();
        b37Chr3Ref = FuncotatorReferenceTestUtils.retrieveB37Chr3Ref();
        hg19Chr3Ref = FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref();
        hg19Chr19Ref = FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref();
    }

    //==================================================================================================================
    // Helper methods to create output files and maybe leave them around to debug the test output.

    private static File getOutputFile(final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType) {
        return getOutputFile("funcotator_tmp_out", outputFormatType.toString().toLowerCase());
    }

    private static File getOutputFile(final String outfileBaseName,
                                      final String outFileExtension) {
        final File outputFile;
        if (!doDebugTests) {
            outputFile = createTempFile(tmpOutDir + File.separator + outfileBaseName, "." + outFileExtension);
        } else {
            outputFile = new File(tmpOutDir, outfileBaseName + "." + outFileExtension);
        }
        return outputFile;
    }

    private static void addManualAnnotationsToArguments(final ArgumentsBuilder arguments) {

        // ================================================================================
        // Annotation Defaults:
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "dbSNP_RS:0");
        arguments.addArgumentWithValueThatIncludesWhitespace(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "dbSNP_Val_Status:No Value");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Center:broad.mit.edu");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "source:WES");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "normal_barcode:normal_sample");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "tumor_barcode:tumor_sample");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "NCBI_Build:37");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Strand:+");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "status:Somatic");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "phase:Phase_I");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "sequencer:Illumina");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Tumor_Validation_Allele1:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Tumor_Validation_Allele2:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Match_Norm_Validation_Allele1:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Match_Norm_Validation_Allele2:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Verification_Status:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Validation_Status:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Validation_Method:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Score:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "BAM_file:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Match_Norm_Seq_Allele1:");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "Match_Norm_Seq_Allele2:");

        // ================================================================================
        // Annotation Overrides:
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME, "Oreganno_Build:BUILDED_GOOD_REAL_BIG");
    }

    //==================================================================================================================

    @DataProvider
    Iterator<Object[]> provideForIntegrationTest() {

        final ArrayList<Object[]> testCases = new ArrayList<>();

        // TODO: Fix this set of tests!  THEY DON'T ALL PASS!  (Issue: https://github.com/broadinstitute/gatk/issues/4344)
        // VCF SNP / MNP / INDEL tests
//        for ( final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values() ) {
//            for ( final SimpleKeyXsvFuncotationFactory.XsvDataKeyType keyType : SimpleKeyXsvFuncotationFactory.XsvDataKeyType.values() ) {
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
//                                FuncotatorTestConstants.PIK3CA_SNP_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.PIK3CA_SNP_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
//                                FuncotatorTestConstants.MUC16_MNP_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.MUC16_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.MUC16_MNP_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//                testCases.add(
//                        new Object[] {
//                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
//                                FuncotatorArgumentDefinitions.ReferenceVersionType.hg19,
//                                FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
//                                FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME + ".vcf",
//                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
//                                keyType,
//                                outputFormatType,
//                                FuncotatorTestConstants.PIK3CA_INDEL_FILE_BASE_NAME + ".oncotatorAnnotated." + outputFormatType.toString().toLowerCase()
//                        }
//                );
//            }
//        }

        // Basic Test Cases:
        for (final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType : FuncotatorArgumentDefinitions.OutputFormatType.values()) {
            for (final SimpleKeyXsvFuncotationFactory.XsvDataKeyType keyType : SimpleKeyXsvFuncotationFactory.XsvDataKeyType.values()) {
                testCases.add(
                        new Object[]{
                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                                hg19Chr3Ref,
                                FuncotatorTestConstants.VARIANT_FILE_HG19_CHR3,
                                FuncotatorTestConstants.PIK3CA_TRANSCRIPT,
                                keyType,
                                outputFormatType,
                                null
                        }
                );
                testCases.add(
                        new Object[]{
                                FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER,
                                FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                                hg19Chr19Ref,
                                FuncotatorTestConstants.VARIANT_FILE_HG19_CHR19,
                                FuncotatorTestConstants.MUC16_TRANSCRIPT,
                                keyType,
                                outputFormatType,
                                null
                        }
                );
            }
        }

        return testCases.iterator();
    }

    @DataProvider
    public Object[][] provideForLargeDataValidationTest() {
        return new Object[][]{
                {
                        "M2_01115161-TA1-filtered.vcf",
                        "Homo_sapiens_assembly19.fasta",
                        true,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19,
                },
                {
                        "C828.TCGA-D3-A2JP-06A-11D-A19A-08.3-filtered.PASS.vcf",
                        "Homo_sapiens_assembly19.fasta",
                        true,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG19
                },
                {
                        "hg38_test_variants.vcf",
                        "Homo_sapiens_assembly38.fasta",
                        false,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38
                },
                {
                        "sample21.trimmed.vcf",
                        "Homo_sapiens_assembly38.fasta",
                        false,
                        FuncotatorTestConstants.REFERENCE_VERSION_HG38
                }
        };
    }

    //==================================================================================================================

    // This test is to make sure we don't create a bunch of temp files anywhere.
    // It will force anyone who changes the outputToTmpDir flag to make it true when they check in this test file.
    @Test(groups = {"funcotatorValidation"})
    public void metaTestEnsureTempDirs() {
        Assert.assertEquals(doDebugTests, false);
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void testFuncotatorWithoutValidatingResults(final String dataSourcesPath,
                                                       final String refVer,
                                                       final String referenceFileName,
                                                       final String variantFileName,
                                                       final String transcriptName,
                                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                                       final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                                       final String expectedOutputFilePath) throws IOException {

        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(variantFileName));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(referenceFileName));
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, dataSourcesPath);
        arguments.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, refVer);
        arguments.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());

        runCommandLine(arguments);
    }

    @Test(enabled = doDebugTests,
            groups = {"funcotatorValidation"},
            dataProvider = "provideForLargeDataValidationTest")
    public void largeDataValidationTest(final String inputVcfName,
                                        final String referencePath,
                                        final boolean allowHg19B37ContigMatches,
                                        final String referenceVersion) throws IOException {

        // Get our main test folder path from our environment:
        final String testFolderInputPath = getFuncotatorLargeDataValidationTestInputPath();

        long startTime = 0, endTime = 0;
        final long overallStartTime = System.nanoTime();

        final String outFileBaseName = inputVcfName + ".funcotator";

        for (final FuncotatorArgumentDefinitions.OutputFormatType outFormat : FuncotatorArgumentDefinitions.OutputFormatType.values()) {

            startTime = System.nanoTime();

            final File outputFile;
            if (outFormat == FuncotatorArgumentDefinitions.OutputFormatType.VCF) {
                outputFile = getOutputFile(outFileBaseName, outFormat.toString().toLowerCase());
            } else {
                outputFile = getOutputFile(outFileBaseName + ".maf", "tsv");
            }

            final ArgumentsBuilder arguments = new ArgumentsBuilder();

            arguments.addArgument("verbosity", "INFO");

            arguments.addVCF(new File(testFolderInputPath + inputVcfName));

            arguments.addReference(new File(testFolderInputPath + referencePath));

            arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, getFuncotatorLargeDataValidationTestInputPath() + LARGE_DATASOURCES_FOLDER);

            arguments.addBooleanArgument(FuncotatorArgumentDefinitions.ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_LONG_NAME, allowHg19B37ContigMatches);

            arguments.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, referenceVersion);
            arguments.addOutput(outputFile);
            arguments.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outFormat.toString());

            // Add our manual annotations to the arguments:
            addManualAnnotationsToArguments(arguments);

            // Run the tool with our args:
            runCommandLine(arguments);

            endTime = System.nanoTime();

            System.out.println("  Elapsed Time (" + outFormat.toString() + "): " + (endTime - startTime) / 1e9 + "s");
        }

        System.out.println("Total Elapsed Time: " + (endTime - overallStartTime) / 1e9 + "s");
    }

    @Test(dataProvider = "provideForIntegrationTest")
    public void exhaustiveArgumentTest(final String dataSourcesPath,
                                       final String refVer,
                                       final String referenceFileName,
                                       final String variantFileName,
                                       final String transcriptName,
                                       final SimpleKeyXsvFuncotationFactory.XsvDataKeyType xsvMatchType,
                                       final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType,
                                       final String expectedOutputFilePath) throws IOException {

        final String outFileName = "funcotator_tmp_out_" + xsvMatchType.toString() + "_" + transcriptName + "." + outputFormatType.toString().toLowerCase();

        final File outputFile = getOutputFile(outFileName.substring(0, outFileName.length() - 4), outFileName.substring(outFileName.length() - 3));

        // Set up our transcript of interest in a file:
        final File transcriptIdFile = getSafeNonExistentFile("TranscriptIdFile.txt");
        try (final BufferedWriter bw = new BufferedWriter(new FileWriter(transcriptIdFile))) {
            bw.write(transcriptName);
        }

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        // Required Args:
        arguments.addVCF(new File(variantFileName));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(referenceFileName));

        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, dataSourcesPath);
        arguments.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, refVer);
        arguments.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());

        // Transcript selection:
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME, TranscriptSelectionMode.BEST_EFFECT.toString());
        arguments.addArgument(FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME, transcriptIdFile.toString());

        // Annotation Defaults and Overrides:
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME, "GARBAGEDAY:SPUMONI");
        arguments.addArgument(FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME, "Oreganno_Build:BUILDED_GOOD_REAL_BIG");

        // Run the beast:
        runCommandLine(arguments);

        // Only test for content-correctness if the output file was specified:
        if (expectedOutputFilePath != null) {
            // Get the expected output file:
            final File expectedOutputFile = new File(expectedOutputFilePath);

            // Make sure that the actual and expected output files are the same:
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedOutputFile, "#");
        }
    }

    /**
     * Test that we can annotate a hg19 datasource when GENCODE is using "chr*" and the datasource is not.
     */
    @Test
    public void testCanAnnotateMixedContigHg19Clinvar() {
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(PIK3CA_VCF_HG19));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(b37Chr3Ref));
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_PIK3CA_DIR);
        arguments.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, FuncotatorTestConstants.REFERENCE_VERSION_HG19);
        arguments.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());
        arguments.addBooleanArgument(FuncotatorArgumentDefinitions.ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_LONG_NAME, true);

        runCommandLine(arguments);

        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        final int NUM_VARIANTS = 21;
        final int NUM_CLINVAR_HITS = 4;
        Assert.assertEquals(variantContexts.size(), NUM_VARIANTS);

        // Look for "MedGen" to know that we have a clinvar hit.
        Assert.assertEquals(variantContexts.stream()
                .filter(vc -> StringUtils.contains(vc.getAttributeAsString("FUNCOTATION", ""), "MedGen"))
                .count(), NUM_CLINVAR_HITS);
    }

    @Test
    public void testCanAnnotateHg38ClinvarAndGencodeV28() {
        // Clinvar datasource did  go through one round of preprocessing to make contig names "1" --> "chr1" (for example).  This is an issue with ClinVar, not GATK.
        final FuncotatorArgumentDefinitions.OutputFormatType outputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File outputFile = getOutputFile(outputFormatType);

        final ArgumentsBuilder arguments = new ArgumentsBuilder();

        arguments.addVCF(new File(PIK3CA_VCF_HG38));
        arguments.addOutput(outputFile);
        arguments.addReference(new File(hg38Chr3Ref));
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_PIK3CA_DIR);
        arguments.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, FuncotatorTestConstants.REFERENCE_VERSION_HG38);
        arguments.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, outputFormatType.toString());

        runCommandLine(arguments);

        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        final int NUM_VARIANTS = 100;
        final int NUM_CLINVAR_HITS = 1; // This was verified with SelectVariants
        Assert.assertEquals(variantContexts.size(), NUM_VARIANTS);

        // Look for "MedGen" to know that we have a clinvar hit.
        Assert.assertEquals(variantContexts.stream()
                .filter(vc -> StringUtils.contains(vc.getAttributeAsString("FUNCOTATION", ""), "MedGen"))
                .count(), NUM_CLINVAR_HITS);
    }

    @DataProvider(name="provideForMafVcfConcordanceProteinChangePik3ca")
    final Object[][] provideForMafVcfConcordanceProteinChangePik3ca() {
        return new Object[][] {
                {PIK3CA_VCF_HG19_SNPS, b37Chr3Ref, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Arrays.asList("Gencode_19_proteinChange"), DS_PIK3CA_DIR, 15},
                {PIK3CA_VCF_HG19_INDELS, b37Chr3Ref, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Arrays.asList("Gencode_19_proteinChange"), DS_PIK3CA_DIR, 57},
                {MUC16_VCF_HG19, hg19Chr19Ref, FuncotatorTestConstants.REFERENCE_VERSION_HG19, Arrays.asList("Gencode_19_proteinChange"), DS_MUC16_DIR, 2057}
        };
    }

    /**
     * Make sure that VCFs and MAFs have exactly the same protein change strings.
     */
    @Test(dataProvider="provideForMafVcfConcordanceProteinChangePik3ca")
    public void testVcfMafConcordanceForProteinChange(final String inputVcf, final String inputRef,
                                                          final String funcotatorRef, final List<String> annotationsToCheck,
                                                          final String datasourceDir,
                                                          final int gtNumVariants) {
        final FuncotatorArgumentDefinitions.OutputFormatType vcfOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File vcfOutputFile = getOutputFile(vcfOutputFormatType);

        final ArgumentsBuilder argumentsVcf = new ArgumentsBuilder();

        argumentsVcf.addVCF(new File(inputVcf));
        argumentsVcf.addOutput(vcfOutputFile);
        argumentsVcf.addReference(new File(inputRef));
        argumentsVcf.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, datasourceDir);
        argumentsVcf.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, funcotatorRef);
        argumentsVcf.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, vcfOutputFormatType.toString());
        argumentsVcf.addBooleanArgument(FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME, false);

        // We need this argument since we are testing on a subset of b37
        argumentsVcf.addBooleanArgument(FuncotatorArgumentDefinitions.ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_OVERRIDE_LONG_NAME, true);
        runCommandLine(argumentsVcf);

        final FuncotatorArgumentDefinitions.OutputFormatType mafOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.MAF;
        final File mafOutputFile = getOutputFile(mafOutputFormatType);

        final ArgumentsBuilder argumentsMaf = new ArgumentsBuilder();

        argumentsMaf.addVCF(new File(inputVcf));
        argumentsMaf.addOutput(mafOutputFile);
        argumentsMaf.addReference(new File(inputRef));
        argumentsMaf.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, datasourceDir);
        argumentsMaf.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, funcotatorRef);
        argumentsMaf.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, mafOutputFormatType.toString());
        argumentsMaf.addBooleanArgument(FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME, false);
        argumentsMaf.addBooleanArgument(FuncotatorArgumentDefinitions.ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_OVERRIDE_LONG_NAME, true);

        runCommandLine(argumentsMaf);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(vcfOutputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);

        Assert.assertTrue(variantContexts.stream().allMatch(v -> v.hasAttribute(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME)));

        final AnnotatedIntervalCollection maf = AnnotatedIntervalCollection.create(mafOutputFile.toPath(), null);
        Assert.assertEquals(maf.getRecords().size(), gtNumVariants);
        Assert.assertTrue(maf.getRecords().stream()
                .anyMatch(v -> !v.getAnnotationValue(MafOutputRendererConstants.FieldName_Variant_Classification).equals("IGR")));
        Assert.assertTrue(maf.getRecords().stream()
                .anyMatch(v -> v.getAnnotationValue(MafOutputRendererConstants.FieldName_Variant_Classification).equals("Missense_Mutation") ||
                        v.getAnnotationValue(MafOutputRendererConstants.FieldName_Variant_Classification).startsWith("Frame_Shift")));

        // Get the protein changes:
        final String[] funcotationKeys = extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());
        for (final String annotationToCheck: annotationsToCheck) {
            final List<String> mafProteinChanges = maf.getRecords().stream().map(v -> v.getAnnotationValue(MafOutputRendererConstants.FieldName_Protein_Change)).collect(Collectors.toList());
            final List<String> vcfProteinChanges = variantContexts.stream().map(v -> v.getAttributeAsString(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME, null))
                    .map(f -> FuncotatorUtils.getFuncotationMapFromVcfFuncotationField(funcotationKeys, f))
                    .map(m -> m.get(annotationToCheck)).collect(Collectors.toList());
            Assert.assertEquals(mafProteinChanges, vcfProteinChanges, "Failed matching " + annotationToCheck);
        }
    }

    @Test
    public void testVcfDatasourceAccountsForAltAlleles() {
        final FuncotatorArgumentDefinitions.OutputFormatType vcfOutputFormatType = FuncotatorArgumentDefinitions.OutputFormatType.VCF;
        final File vcfOutputFile = getOutputFile(vcfOutputFormatType);

        final ArgumentsBuilder argumentsVcf = new ArgumentsBuilder();

        argumentsVcf.addVCF(new File(PIK3CA_VCF_HG19_ALTS));
        argumentsVcf.addOutput(vcfOutputFile);
        argumentsVcf.addReference(new File(FuncotatorTestConstants.HG19_3_REFERENCE_FILE_NAME));
        argumentsVcf.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_PIK3CA_DIR);
        argumentsVcf.addArgument(FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME, FuncotatorTestConstants.REFERENCE_VERSION_HG19);
        argumentsVcf.addArgument(FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME, vcfOutputFormatType.toString());
        argumentsVcf.addBooleanArgument(FuncotatorArgumentDefinitions.REMOVE_FILTERED_VARIANTS_LONG_NAME, false);

        // We need this argument since we are testing on a subset of b37
        argumentsVcf.addBooleanArgument(FuncotatorArgumentDefinitions.ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_OVERRIDE_LONG_NAME, true);
        runCommandLine(argumentsVcf);

        final Pair<VCFHeader, List<VariantContext>> vcfInfo = VariantContextTestUtils.readEntireVCFIntoMemory(vcfOutputFile.getAbsolutePath());
        final List<VariantContext> variantContexts = vcfInfo.getRight();
        final VCFHeader vcfHeader = vcfInfo.getLeft();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        final String[] funcotationKeys = extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());

        // The first variant context should have clinvar annotations, since it hit on the alt allele.  None of the rest.
        // This test assumes that each test variant context has only one alt allele.
        final String funcotationInfoFieldWithClinVarHit = variantContexts.get(0).getAttributeAsString(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME, null);
        Assert.assertEquals(FuncotatorUtils.getFuncotationMapFromVcfFuncotationField(funcotationKeys, funcotationInfoFieldWithClinVarHit).get("dummy_ClinVar_VCF_CLNDISDB"),
            FuncotatorUtils.sanitizeFuncotationForVcf("MedGen:C0027672,SNOMED_CT:699346009"));

        // The rest should not have any clinvar hits.
        final List<String> clinvarAnnotations = new ArrayList<>();
        final List<String> clinvarAnnotationsTruth = new ArrayList<>();
        for (int i = 1; i < variantContexts.size(); i++) {
            final String infoField = variantContexts.get(i).getAttributeAsString(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME, null);
            clinvarAnnotations.add(FuncotatorUtils.getFuncotationMapFromVcfFuncotationField(funcotationKeys, infoField).get("dummy_ClinVar_VCF_CLNDISDB"));
            clinvarAnnotationsTruth.add(FuncotatorUtils.sanitizeFuncotationForVcf(""));
        }
        Assert.assertEquals(clinvarAnnotations, clinvarAnnotationsTruth);
    }
}
