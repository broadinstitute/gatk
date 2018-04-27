package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledModeledSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledModeledSegment;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;
import java.util.List;
import java.io.File;
import java.util.Scanner;

public class CallModeledSegmentsIntegrationTest extends CommandLineProgramTest {
    // Simulated samples
    private static final File SIMULATED_DATA_DIR = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/");
    private static final File INPUT_SIMULATED_NORMAL_DATA = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/normal_sample_data.seg");
    private static final File INPUT_SIMULATED_NORMAL_TRUTH = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/normal_sample_truth.seg");
    private static final File INPUT_SIMULATED_NORMAL_NOISY_DATA = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/normal_noisy_sample_data.seg");
    private static final File INPUT_SIMULATED_NORMAL_NOISY_TRUTH = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/normal_noisy_sample_truth.seg");
    private static final File INPUT_SIMULATED_CANCER_40P_PURITY_DATA = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/cancer_0p40_purity_data.seg");
    private static final File INPUT_SIMULATED_CANCER_40P_PURITY_TRUTH = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/cancer_0p40_purity_truth.seg");
    private static final File INPUT_SIMULATED_CANCER_60P_PURITY_DATA = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/cancer_0p60_purity_data.seg");
    private static final File INPUT_SIMULATED_CANCER_60P_PURITY_TRUTH = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/cancer_0p60_purity_truth.seg");
    private static final File INPUT_SIMULATED_CANCER_100P_PURITY_DATA = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/cancer_1p00_purity_data.seg");
    private static final File INPUT_SIMULATED_CANCER_100P_PURITY_TRUTH = new File("/Users/mkanaszn/Broad_Institute/Code/gatk_ssh/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/modeled-segments-caller-sim-data/cancer_1p00_purity_truth.seg");

    // Real samples
    // private static final File INPUT_NORMAL = new File(
    //         "/Volumes/dsde_working/slee/archived/2017/ModelSegments-pipeline-test-cfc/run/cromwell-executions/CNVSomaticPairsWorkflow/99e77cf5-1aaf-4517-8926-8e0cb38c3222/call-CNVSomaticPairWorkflow/shard-12/CNVSomaticPairWorkflow/3adbf26d-78e0-4da0-96a7-ac48bf5b8648/call-ModelSegmentsNormal/execution/TCGA-44-6147-10A-01D-1753-08.modelFinal.seg");
    // private static final File INPUT_TUMOR = new File(
    //         "/Volumes/dsde_working/slee/archived/2017/ModelSegments-pipeline-test-cfc/run/cromwell-executions/CNVSomaticPairsWorkflow/99e77cf5-1aaf-4517-8926-8e0cb38c3222/call-CNVSomaticPairWorkflow/shard-14/CNVSomaticPairWorkflow/d612e80d-bb5f-4db2-bf86-3b3ce296dafd/call-ModelSegmentsTumor/execution/TCGA-44-6774-01A-21D-1855-08.modelFinal.seg");
    // private static final File INPUT_CELL_LINE_50_PER_CENT_PURITY = new File(
    //         "/Users/mkanaszn/Broad_Institute/Code/Tumor_segmentation/Testing_Lee_Lischtenstein_data/SM-74P3M.modelFinal.seg");

    @Test
    public void testTmp() {
        final File FigOutputDir = createTempDir(SIMULATED_DATA_DIR.getAbsolutePath() + "test/figs");
        final File FileOutputDir = createTempDir(SIMULATED_DATA_DIR.getAbsolutePath() + "test/calls");
        final File LogOutputDir =  createTempDir(SIMULATED_DATA_DIR.getAbsolutePath() + "test/logs");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, INPUT_SIMULATED_NORMAL_DATA.getAbsolutePath())
                //.addArgument(CallModeledSegments.OUTPUT_IMAGE_DIR_LONG_NAME, FigOutputDir.getAbsolutePath())
                //.addArgument(CallModeledSegments.OUTPUT_CALLS_DIR_LONG_NAME, FileOutputDir.getAbsolutePath())
                //.addArgument(CallModeledSegments.LOG_FILENAME_PREFIX_LONG_NAME, LogOutputDir.getAbsolutePath())
                .addArgument(CallModeledSegments.OUTPUT_IMAGE_DIR_LONG_NAME, "/Users/mkanaszn/Desktop/tmp/figs")
                .addArgument(CallModeledSegments.OUTPUT_CALLS_DIR_LONG_NAME, "/Users/mkanaszn/Desktop/tmp/calls")
                .addArgument(CallModeledSegments.LOG_FILENAME_PREFIX_LONG_NAME, "/Users/mkanaszn/Desktop/tmp/logging/test.log")
                .addArgument(CallModeledSegments.OUTPUT_IMAGE_PREFIX_LONG_NAME, "testFig")
                .addArgument(CallModeledSegments.OUTPUT_CALLS_PREFIX_LONG_NAME, "testCalls")
                .addArgument(CallModeledSegments.LOAD_COPY_RATIO_LONG_NAME, "true")
                .addArgument(CallModeledSegments.LOAD_ALLELE_FRACTION_LONG_NAME, "true")
                .addArgument(CallModeledSegments.INTERACTIVE_RUN_LONG_NAME, "true");
        List<String> args = argsBuilder.getArgsList();
        for (String ar:args) {
            System.out.println(ar);
        }
        runCommandLine(argsBuilder);
        // assertOutputFiles(outputDir);
        File outputFile = new File("/Users/mkanaszn/Desktop/testing/results/normal_data.called.seg");
        File truthFile = new File("/Users/mkanaszn/Desktop/testing/normal_truth.seg");
        compareCallFiles(outputFile, truthFile);


        File calls_file = new File("/Users/mkanaszn/Desktop/tmp/calls/testCalls.called.seg");
        try {
            Scanner sc = new Scanner(calls_file);
            while(sc.hasNextLine()) {
                System.out.println(sc.nextLine());
            }
        } catch(java.io.FileNotFoundException e) {
            System.out.println("File not found!");
        }
    }

    @Test
    public void testNormal() {
        final File FigOutputDir = createTempDir("test/figs");
        final File FileOutputDir = createTempDir("test/calls");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(INPUT_SIMULATED_NORMAL_DATA)
                .addArgument(CallModeledSegments.OUTPUT_IMAGE_DIR_LONG_NAME, FigOutputDir.getAbsolutePath())
                .addArgument(CallModeledSegments.OUTPUT_CALLS_DIR_LONG_NAME, FileOutputDir.getAbsolutePath())
                .addArgument(CallModeledSegments.LOAD_COPY_RATIO_LONG_NAME, "true")
                .addArgument(CallModeledSegments.LOAD_ALLELE_FRACTION_LONG_NAME, "true");
        runCommandLine(argsBuilder);
        //assertOutputFiles(outputDir);
    }


    private static void assertOutputFiles(final File inputFile,
                                          final File figOutputDir,
                                          final File callsOutputFile) {
        Assert.assertTrue(inputFile.isFile());
        Assert.assertTrue(figOutputDir.isDirectory());
        Assert.assertTrue(callsOutputFile.isDirectory());

    }

    private static boolean compareCallFiles(final File outputFile, final File truthFile) {
        CalledModeledSegmentCollection outputData = new CalledModeledSegmentCollection(outputFile);
        CalledModeledSegmentCollection truthData = new CalledModeledSegmentCollection(truthFile);
        List<CalledModeledSegment> outputDataRecords = outputData.getRecords();
        for (int i=0; i<10; i++) {
            System.out.println(outputDataRecords.get(i).getCallNormal());
        }
        //for (CalledModeledSegment o : outputData) {
        //    System.out.println(o.getCallNormal());
        //}


        return true;
    }
}


//    private static void assertOutputFiles(final File outputDir,
//                                          final String outputPrefix,
//                                          final boolean isAllelicCountsPresent,
//                                          final boolean isNormalAllelicCountsPresent) {
//        Assert.assertTrue(!(!isAllelicCountsPresent && isNormalAllelicCountsPresent));
//        for (final String fileTag : Arrays.asList(ModelSegments.BEGIN_FIT_FILE_TAG, ModelSegments.FINAL_FIT_FILE_TAG)) {
//            final ModeledSegmentCollection modeledSegments =
//                    new ModeledSegmentCollection(new File(outputDir, outputPrefix + fileTag + ModelSegments.SEGMENTS_FILE_SUFFIX));
//            Assert.assertEquals(EXPECTED_METADATA, modeledSegments.getMetadata());
//
//            final ParameterDecileCollection<CopyRatioParameter> copyRatioParameters = new ParameterDecileCollection<>(new File(outputDir,
//                    outputPrefix + fileTag + ModelSegments.COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX), CopyRatioParameter.class);
//            Assert.assertEquals(EXPECTED_METADATA.getSampleName(), copyRatioParameters.getMetadata().getSampleName());
//
//            final ParameterDecileCollection<AlleleFractionParameter> alleleFractionParameters = new ParameterDecileCollection<>(new File(outputDir,
//                    outputPrefix + fileTag + ModelSegments.ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX), AlleleFractionParameter.class);
//            Assert.assertEquals(EXPECTED_METADATA.getSampleName(), alleleFractionParameters.getMetadata().getSampleName());
//        }
//
//        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(new File(outputDir,
//                outputPrefix + ModelSegments.COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX));
//        Assert.assertEquals(EXPECTED_METADATA, copyRatioSegments.getMetadata());
//
//        AllelicCountCollection hetAllelicCounts = null;
//        AllelicCountCollection hetNormalAllelicCounts = null;
//        if (isAllelicCountsPresent) {
//            hetAllelicCounts = new AllelicCountCollection(new File(outputDir,
//                    outputPrefix + ModelSegments.HET_ALLELIC_COUNTS_FILE_SUFFIX));
//            Assert.assertEquals(EXPECTED_METADATA, hetAllelicCounts.getMetadata());
//        }
//        if (isNormalAllelicCountsPresent) { //if this is true, case sample allelic counts will be present
//            hetNormalAllelicCounts = new AllelicCountCollection(new File(outputDir,
//                    outputPrefix + ModelSegments.NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX));
//            Assert.assertNotEquals(EXPECTED_METADATA, hetNormalAllelicCounts.getMetadata());    //sample names should differ
//            Assert.assertTrue(CopyNumberArgumentValidationUtils.isSameDictionary(               //sequence dictionary should be the same
//                    EXPECTED_METADATA.getSequenceDictionary(), hetNormalAllelicCounts.getMetadata().getSequenceDictionary()));
//            Assert.assertEquals(hetAllelicCounts.getIntervals(), hetNormalAllelicCounts.getIntervals());
//        }
//    }