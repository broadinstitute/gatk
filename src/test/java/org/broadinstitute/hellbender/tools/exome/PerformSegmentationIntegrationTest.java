package org.broadinstitute.hellbender.tools.exome;

import com.google.common.collect.Lists;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;
import org.broadinstitute.hellbender.utils.segmenter.SegmenterUnitTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class PerformSegmentationIntegrationTest extends CommandLineProgramTest{

    private static final String inputTestDir = "src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/";
    private static final String outputTestDir = "src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/";

    @DataProvider(name="inputFileData")
    public Object[][] inputFileData() {
        return new Object[][] {
                new Object[] { new File(inputTestDir, "HCC1143_reduced.tsv"), new File(outputTestDir, "HCC1143_reduced_result.seg"), createTempFile("gatkcnv.HCC1143", ".seg"), "HCC1143"},
                new Object[] { new File(inputTestDir, "HCC1143_short.tsv"), new File(outputTestDir, "HCC1143_short_result.seg"), createTempFile("gatkcnv.HCC1143", ".seg"), "HCC1143"},
                new Object[] { new File(inputTestDir, "Simple.tsv"), new File(outputTestDir, "Simple_result.seg"), createTempFile("gatkcnv.HCC1143", ".seg"), "Simple"},
        };
    }

    @Test(dataProvider = "inputFileData")
    public void testUnLoggedCommandLine(final File INPUT_FILE, final File EXPECTED, final File output, String sampleName) throws IOException {
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output.getAbsolutePath(),
        };
        runCommandLine(arguments);
        SegmenterUnitTest.assertEqualSegments(output, EXPECTED);
    }

    @Test()
    public void testUnLoggedCommandLine() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "HCC1143_reduced_log.tsv");
        final File EXPECTED = new File(outputTestDir, "HCC1143_reduced_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), true);
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
        };
        runCommandLine(arguments);
        SegmenterUnitTest.assertEqualSegments(output, EXPECTED);
    }

    @Test()
    public void testUnLoggedCommandLineWithWeights() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "HCC1143_reduced_log.tsv");
        final File EXPECTED = new File(outputTestDir, "HCC1143_reduced_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final File tmpWeightsFile = IOUtils.createTempFile("integration-weight-file", ".txt");
        final double [] weights = new double[7677];
        Arrays.fill(weights, 1.0);
        ParamUtils.writeValuesToFile(weights, tmpWeightsFile);
        final String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), true);
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
                "-" + PerformSegmentation.TARGET_WEIGHT_FILE_SHORT_NAME, tmpWeightsFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        SegmenterUnitTest.assertEqualSegments(output, EXPECTED);
    }

    @Test(dataProvider = "parameterTests")
    public void testAlternateParametersActuallyChangeData(String [] newArguments) {

        // This test simply tests that if we change the value of the parameter that the output is altered.

        final File INPUT_FILE = new File(inputTestDir, "HCC1143_reduced_log.tsv");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final File outputNewParam = createTempFile("gatkcnv.HCC1143.sdundo", ".seg");
        final File EXPECTED = new File(outputTestDir, "HCC1143_reduced_result.seg");
        final File tmpWeightsFile = IOUtils.createTempFile("integration-weight-file", ".txt");
        final double [] weights = new double[7677];
        Arrays.fill(weights, 1.0);
        ParamUtils.writeValuesToFile(weights, tmpWeightsFile);
        final String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), true);
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
                "-" + PerformSegmentation.TARGET_WEIGHT_FILE_SHORT_NAME, tmpWeightsFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output.getAbsolutePath(),
        };
        runCommandLine(arguments);
        SegmenterUnitTest.assertEqualSegments(output, EXPECTED);

        final List<String> fullNewArgumentsAsList = Lists.newArrayList(arguments);
        // Change the output file
        fullNewArgumentsAsList.remove("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        fullNewArgumentsAsList.remove(output.getAbsolutePath());
        fullNewArgumentsAsList.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        fullNewArgumentsAsList.add(outputNewParam.getAbsolutePath());
        fullNewArgumentsAsList.addAll(Arrays.asList(newArguments));

        runCommandLine(fullNewArgumentsAsList.toArray(new String[fullNewArgumentsAsList.size()]));
        SegmenterUnitTest.assertNotEqualSegments(outputNewParam, output);
    }

    @Test
    public void testAlternateSDUndoParametersActuallyChangeData() {

        // This test simply tests that if we change the value of the parameter that the output is altered.

        final File INPUT_FILE = new File(inputTestDir, "HCC1143_reduced_log.tsv");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final File outputNewParam = createTempFile("gatkcnv.HCC1143.sdundo", ".seg");

        final File tmpWeightsFile = IOUtils.createTempFile("integration-weight-file", ".txt");
        final double [] weights = new double[7677];
        Arrays.fill(weights, 1.0);
        ParamUtils.writeValuesToFile(weights, tmpWeightsFile);
        final String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), true);
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
                "-" + PerformSegmentation.TARGET_WEIGHT_FILE_SHORT_NAME, tmpWeightsFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output.getAbsolutePath(),
                "-" + PerformSegmentation.UNDOSPLITS_SHORT_NAME, StringUtils.upperCase(RCBSSegmenter.UndoSplits.SDUNDO.toString())
        };
        runCommandLine(arguments);

        final String[] newArguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
                "-" + PerformSegmentation.TARGET_WEIGHT_FILE_SHORT_NAME, tmpWeightsFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputNewParam.getAbsolutePath(),
                "-" + PerformSegmentation.UNDOSPLITS_SHORT_NAME, StringUtils.upperCase(RCBSSegmenter.UndoSplits.SDUNDO.toString()),
                "-" + PerformSegmentation.UNDOSD_SHORT_NAME, String.valueOf(2)
        };
        runCommandLine(newArguments);

        SegmenterUnitTest.assertNotEqualSegments(outputNewParam, output);
    }

    @Test
    public void testAlternateSDUndoPruneParametersActuallyChangeData() {

        // This test simply tests that if we change the value of the parameter that the output is altered.

        final File INPUT_FILE = new File(inputTestDir, "HCC1143_reduced_log.tsv");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final File outputNewParam = createTempFile("gatkcnv.HCC1143.sdundo", ".seg");

        final File tmpWeightsFile = IOUtils.createTempFile("integration-weight-file", ".txt");
        final double [] weights = new double[7677];
        Arrays.fill(weights, 1.0);
        ParamUtils.writeValuesToFile(weights, tmpWeightsFile);
        final String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), true);
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
                "-" + PerformSegmentation.TARGET_WEIGHT_FILE_SHORT_NAME, tmpWeightsFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, output.getAbsolutePath(),
                "-" + PerformSegmentation.UNDOSPLITS_SHORT_NAME, StringUtils.upperCase(RCBSSegmenter.UndoSplits.PRUNE.toString())
        };
        runCommandLine(arguments);

        final String[] newArguments = {
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
                "-" + PerformSegmentation.TARGET_WEIGHT_FILE_SHORT_NAME, tmpWeightsFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputNewParam.getAbsolutePath(),
                "-" + PerformSegmentation.UNDOSPLITS_SHORT_NAME, StringUtils.upperCase(RCBSSegmenter.UndoSplits.PRUNE.toString()),
                "-" + PerformSegmentation.UNDOPRUNE_SHORT_NAME, String.valueOf(0.1)
        };
        runCommandLine(newArguments);

        SegmenterUnitTest.assertNotEqualSegments(outputNewParam, output);
    }

    @DataProvider(name="parameterTests")
    public Object[][] someTargetsHDF5PoNCreationData() {
        return new Object[][] {
                // Parameter "trim" had to be manually verified (that it was being passed in), since no output values
                //  actually changed, causing this test to fail.
                {
                    new String[]{"-" + PerformSegmentation.UNDOSPLITS_SHORT_NAME, StringUtils.upperCase(RCBSSegmenter.UndoSplits.SDUNDO.toString())}
                }, {
                    new String[]{"-" + PerformSegmentation.UNDOSPLITS_SHORT_NAME, StringUtils.upperCase(RCBSSegmenter.UndoSplits.PRUNE.toString())}
                }, {
                    new String[]{"-" + PerformSegmentation.ALPHA_SHORT_NAME, String.valueOf(0.001)}
                }, {
                    new String[]{"-" + PerformSegmentation.PMETHOD_SHORT_NAME, StringUtils.upperCase(RCBSSegmenter.PMethod.PERM.toString())}
                }, {
                    new String[]{"-" + PerformSegmentation.NMIN_SHORT_NAME, String.valueOf(800)}
                }, {
                    new String[]{"-" + PerformSegmentation.KMAX_SHORT_NAME, String.valueOf(10)}
                }, {
                    new String[]{"-" + PerformSegmentation.MINWIDTH_SHORT_NAME, String.valueOf(5)}
                }, {
                    new String[]{"-" + PerformSegmentation.ETA_SHORT_NAME, String.valueOf(0.1)}
            }
        };
    }
}