package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledModeledSegmentCollection;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;
import java.util.List;
import java.io.File;
import java.util.Scanner;

public class CNVCallerIntegrationTest extends CommandLineProgramTest {
    private static final File INPUT_NORMAL = new File(
            "/Volumes/dsde_working/slee/archived/2017/ModelSegments-pipeline-test-cfc/run/cromwell-executions/CNVSomaticPairsWorkflow/99e77cf5-1aaf-4517-8926-8e0cb38c3222/call-CNVSomaticPairWorkflow/shard-12/CNVSomaticPairWorkflow/3adbf26d-78e0-4da0-96a7-ac48bf5b8648/call-ModelSegmentsNormal/execution/TCGA-44-6147-10A-01D-1753-08.modelFinal.seg");
    private static final File INPUT_TUMOR = new File(
            "/Volumes/dsde_working/slee/archived/2017/ModelSegments-pipeline-test-cfc/run/cromwell-executions/CNVSomaticPairsWorkflow/99e77cf5-1aaf-4517-8926-8e0cb38c3222/call-CNVSomaticPairWorkflow/shard-14/CNVSomaticPairWorkflow/d612e80d-bb5f-4db2-bf86-3b3ce296dafd/call-ModelSegmentsTumor/execution/TCGA-44-6774-01A-21D-1855-08.modelFinal.seg");
    private static final File INPUT_CELL_LINE_50_PER_CENT_PURITY = new File(
            "/Users/mkanaszn/Broad_Institute/Code/Tumor_segmentation/Testing_Lee_Lischtenstein_data/SM-74P3M.modelFinal.seg");

    @Test
    public void testTmp() {
        final File FigOutputDir = createTempDir("/Users/mkanaszn/Desktop/test/figs");
        final File FileOutputDir = createTempDir("/Users/mkanaszn/Desktop/test/calls");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, INPUT_NORMAL.getAbsolutePath())
                //.addArgument(CNVCaller.OUTPUT_IMAGE_DIR_LONG_NAME, FigOutputDir.getAbsolutePath())
                //.addArgument(CNVCaller.OUTPUT_CALLS_DIR_LONG_NAME, FileOutputDir.getAbsolutePath())
                .addArgument(CNVCaller.OUTPUT_IMAGE_DIR_LONG_NAME, "/Users/mkanaszn/Desktop/tmp/figs")
                .addArgument(CNVCaller.OUTPUT_CALLS_DIR_LONG_NAME, "/Users/mkanaszn/Desktop/tmp/calls")
                .addArgument(CNVCaller.OUTPUT_IMAGE_PREFIX_LONG_NAME, "testFig")
                .addArgument(CNVCaller.OUTPUT_CALLS_PREFIX_LONG_NAME, "testCalls")
                .addArgument(CNVCaller.LOAD_COPY_RATIO_LONG_NAME, "true")
                .addArgument(CNVCaller.LOAD_ALLELE_FRACTION_LONG_NAME, "true")
                .addArgument(CNVCaller.INTERACTIVE_RUN_LONG_NAME, "true")
                .addArgument(CNVCaller.LOG_FILENAME_PREFIX_LONG_NAME, "/Users/mkanaszn/Desktop/tmp/logging/test.log");
        List<String> args = argsBuilder.getArgsList();
        for (String ar:args) {
            System.out.println(ar);
        }
        runCommandLine(argsBuilder);
        //assertOutputFiles(outputDir);


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
                .addInput(INPUT_NORMAL)
                .addArgument(CNVCaller.OUTPUT_IMAGE_DIR_LONG_NAME, FigOutputDir.getAbsolutePath())
                .addArgument(CNVCaller.OUTPUT_CALLS_DIR_LONG_NAME, FileOutputDir.getAbsolutePath())
                .addArgument(CNVCaller.LOAD_COPY_RATIO_LONG_NAME, "true")
                .addArgument(CNVCaller.LOAD_ALLELE_FRACTION_LONG_NAME, "true");
        runCommandLine(argsBuilder);
        //assertOutputFiles(outputDir);
    }

    @Test
    public void testTumor() {
        final File FigOutputDir = createTempDir("test/figs");
        final File FileOutputDir = createTempDir("test/calls");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(INPUT_TUMOR)
                .addArgument(CNVCaller.OUTPUT_IMAGE_DIR_LONG_NAME, FigOutputDir.getAbsolutePath())
                .addArgument(CNVCaller.OUTPUT_CALLS_DIR_LONG_NAME, FileOutputDir.getAbsolutePath())
                .addArgument(CNVCaller.LOAD_COPY_RATIO_LONG_NAME, "true")
                .addArgument(CNVCaller.LOAD_ALLELE_FRACTION_LONG_NAME, "true");
        runCommandLine(argsBuilder);
        //assertOutputFiles(outputDir);
    }

    @Test
    public void testCellLine() {
        final File FigOutputDir = createTempDir("test/figs");
        final File FileOutputDir = createTempDir("test/calls");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(INPUT_CELL_LINE_50_PER_CENT_PURITY)
                .addArgument(CNVCaller.OUTPUT_IMAGE_DIR_LONG_NAME, FigOutputDir.getAbsolutePath())
                .addArgument(CNVCaller.OUTPUT_CALLS_DIR_LONG_NAME, FileOutputDir.getAbsolutePath())
                .addArgument(CNVCaller.LOAD_COPY_RATIO_LONG_NAME, "true")
                .addArgument(CNVCaller.LOAD_ALLELE_FRACTION_LONG_NAME, "true");
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