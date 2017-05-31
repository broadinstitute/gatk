package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;


public class SegmenterUnitTest extends BaseTest {

    private static final String inputTestDir = "src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/";
    private static final String outputTestDir = "src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/";

    @Test
    public void testHCC1143Reduced() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "HCC1143_reduced.tsv");
        final File EXPECTED = new File(outputTestDir, "HCC1143_reduced_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testHCC1143Short() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "HCC1143_short.tsv");
        final File EXPECTED = new File(outputTestDir, "HCC1143_short_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testSimple() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "Simple.tsv");
        final File EXPECTED = new File(outputTestDir, "Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "Simple";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testSimpleWithWeights() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "Simple.tsv");
        final File EXPECTED = new File(outputTestDir, "Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "Simple";
        final double [] weights = new double[20];
        Arrays.fill(weights, 1.0);
        final File weightsTmpFile = IOUtils.createTempFile("weights-simple", ".txt");
        ParamUtils.writeValuesToFile(weights, weightsTmpFile);
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false, weightsTmpFile);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testSimpleWithWeightsBig() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "Simple.tsv");
        final File EXPECTED = new File(outputTestDir, "Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "Simple";
        final double [] weights = new double[20];
        Arrays.fill(weights, 10.0);
        weights[10] = 200;
        final File weightsTmpFile = IOUtils.createTempFile("weights-simple", ".txt");
        ParamUtils.writeValuesToFile(weights, weightsTmpFile);
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false, weightsTmpFile);
        assertEqualSegments(output,EXPECTED);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSimpleWithWeightsPosInf() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "Simple.tsv");
        final File EXPECTED = new File(outputTestDir, "Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "Simple";
        final double [] weights = new double[20];
        Arrays.fill(weights, 1.0);
        weights[10] = Double.POSITIVE_INFINITY;
        final File weightsTmpFile = IOUtils.createTempFile("weights-simple", ".txt");
        ParamUtils.writeValuesToFile(weights, weightsTmpFile);
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false, weightsTmpFile);
        assertEqualSegments(output,EXPECTED);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSimpleWithWeightsNegInf() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "Simple.tsv");
        final File EXPECTED = new File(outputTestDir, "Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "Simple";
        final double [] weights = new double[20];
        Arrays.fill(weights, 1.0);
        weights[10] = Double.NEGATIVE_INFINITY;
        final File weightsTmpFile = IOUtils.createTempFile("weights-simple", ".txt");
        ParamUtils.writeValuesToFile(weights, weightsTmpFile);
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), true, weightsTmpFile);
        assertEqualSegments(output,EXPECTED);
    }


    @Test(expectedExceptions = GATKException.class)
    public void testSimpleWithWeightsNegInfNotLogged() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "Simple.tsv");
        final File EXPECTED = new File(outputTestDir, "Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "Simple";
        final double [] weights = new double[20];
        Arrays.fill(weights, 1.0);
        weights[10] = Double.NEGATIVE_INFINITY;
        final File weightsTmpFile = IOUtils.createTempFile("weights-simple", ".txt");
        ParamUtils.writeValuesToFile(weights, weightsTmpFile);
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false, weightsTmpFile);
        assertEqualSegments(output,EXPECTED);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSimpleWithWeightsNan() throws IOException {
        final File INPUT_FILE = new File(inputTestDir, "Simple.tsv");
        final File EXPECTED = new File(outputTestDir, "Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        final String sampleName = "Simple";
        final double [] weights = new double[20];
        Arrays.fill(weights, 1.0);
        weights[10] = Double.NaN;
        final File weightsTmpFile = IOUtils.createTempFile("weights-simple", ".txt");
        ParamUtils.writeValuesToFile(weights, weightsTmpFile);
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false, weightsTmpFile);
        assertEqualSegments(output,EXPECTED);
    }

    /**
     * Compares the content of two segmenter output files.
     * @param actualOutput the actual segmenter output containing file.
     * @param expectedOutput the expected segmenter output containing file.
     * @throws NullPointerException if either {@code actualOutput} or {@code expectedOutput} is {@code null}.
     * @throws AssertionError if there are differences in the segments between both files.
     */
    public static void assertEqualSegments(final File actualOutput, final File expectedOutput) {
        final List<ModeledSegment> actualSegments = SegmentUtils.readModeledSegmentsFromSegmentFile(actualOutput);
        final List<ModeledSegment> expectedSegments = SegmentUtils.readModeledSegmentsFromSegmentFile(expectedOutput);
        Assert.assertEquals(actualSegments, expectedSegments);
    }

    /**
     * Compares the content of two segmenter output files and makes sure that they are NOT the same..
     * @param left the actual segmenter output containing file.
     * @param right the expected segmenter output containing file.
     * @throws NullPointerException if either {@code left} or {@code right} is {@code null}.
     * @throws AssertionError if there are no differences in the segments between both files.
     */
    public static void assertNotEqualSegments(final File left, final File right) {
        final List<ModeledSegment> actualSegments = SegmentUtils.readModeledSegmentsFromSegmentFile(left);
        final List<ModeledSegment> expectedSegments = SegmentUtils.readModeledSegmentsFromSegmentFile(right);
        Assert.assertNotEquals(actualSegments, expectedSegments);
    }
}