package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;
import org.broadinstitute.hellbender.utils.segmenter.SegmenterUnitTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class PerformSegmentationIntegrationTest extends CommandLineProgramTest{

    @DataProvider(name="inputFileData")
    public Object[][] inputFileData() {
        return new Object[][] {
                new Object[] { new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/HCC1143_reduced.tsv"), new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/HCC1143_reduced_result.seg"), createTempFile("recapseg.HCC1143", ".seg"), "HCC1143"},
                new Object[] { new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/HCC1143_short.tsv"), new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/HCC1143_short_result.seg"), createTempFile("recapseg.HCC1143", ".seg"), "HCC1143"},
                new Object[] { new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/Simple.tsv"), new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/Simple_result.seg"), createTempFile("recapseg.HCC1143", ".seg"), "Simple"},
        };
    }

    @Test(dataProvider = "inputFileData")
    public void testUnLoggedCommandLine(final File INPUT_FILE, final File EXPECTED, final File output, String sampleName) throws IOException {
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        final String[] arguments = {
                "-" + PerformSegmentation.SAMPLE_NAME_SHORT_NAME, sampleName,
                "-" + PerformSegmentation.TARGETS_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + PerformSegmentation.SEGMENT_FILE_SHORT_NAME, output.getAbsolutePath(),
        };
        runCommandLine(arguments);
        SegmenterUnitTest.assertEqualSegments(output, EXPECTED);
    }

    @Test()
    public void testUnLoggedCommandLine() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/HCC1143_reduced_log.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/HCC1143_reduced_result.seg");
        final File output = createTempFile("recapseg.HCC1143", ".seg");
        final String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), true);
        final String[] arguments = {
                "-" + PerformSegmentation.SAMPLE_NAME_SHORT_NAME, sampleName,
                "-" + PerformSegmentation.TARGETS_FILE_SHORT_NAME, INPUT_FILE.getAbsolutePath(),
                "-" + PerformSegmentation.SEGMENT_FILE_SHORT_NAME, output.getAbsolutePath(),
                "-" + PerformSegmentation.LOG2_SHORT_NAME,
        };
        runCommandLine(arguments);
        SegmenterUnitTest.assertEqualSegments(output, EXPECTED);
    }
}