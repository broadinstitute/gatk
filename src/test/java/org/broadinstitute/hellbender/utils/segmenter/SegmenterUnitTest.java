package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


public class SegmenterUnitTest extends BaseTest {

    @Test
    public void testHCC1143Reduced() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/HCC1143_reduced.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/HCC1143_reduced_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testHCC1143Short() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/HCC1143_short.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/HCC1143_short_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testSimple() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "Simple";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testSimpleWithWeights() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "Simple";
        final double [] weights = new double[20];
        Arrays.fill(weights, 1.0);
        final File weightsTmpFile = IOUtils.createTempFile("weights-simple", ".txt");
        ParamUtils.writeValuesToFile(weights, weightsTmpFile);
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false, weightsTmpFile);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testSimpleWithWeightsBig() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "Simple";
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
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "Simple";
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
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "Simple";
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
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "Simple";
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
        final File INPUT_FILE = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/org/broadinstitute/hellbender/utils/segmenter/output/Simple_result.seg");
        final File output = createTempFile("gatkcnv.HCC1143", ".seg");
        String sampleName = "Simple";
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
     * @throws AssertionError if there are significant between both files.
     */
    public static void assertEqualSegments(final File actualOutput, final File expectedOutput) {
        try (final SegmentReader actual = new SegmentReader(actualOutput);
             final SegmentReader expected = new SegmentReader(expectedOutput)) {
            final List<SegmentMean> actualSegmentMeans = actual.stream().collect(Collectors.toList());
            final List<SegmentMean> expectedSegmentMeans = expected.stream().collect(Collectors.toList());
            Assert.assertEquals(actualSegmentMeans, expectedSegmentMeans);
        } catch (final IOException ioe) {
            Assert.fail("Unit test threw exception trying to read segments.", ioe);
        }
    }

    /**
     * Compares the content of two segmenter output files and makes sure that they are NOT the same..
     * @param left the actual segmenter output containing file.
     * @param right the expected segmenter output containing file.
     * @throws NullPointerException if either {@code left} or {@code right} is {@code null}.
     * @throws AssertionError if the files are the same between both files.
     */
    public static void assertNotEqualSegments(final File left, final File right) {
        try (final SegmentReader actual = new SegmentReader(left);
             final SegmentReader expected = new SegmentReader(right)) {
            final List<SegmentMean> actualSegmentMeans = actual.stream().collect(Collectors.toList());
            final List<SegmentMean> expectedSegmentMeans = expected.stream().collect(Collectors.toList());
            Assert.assertNotEquals(actualSegmentMeans, expectedSegmentMeans, "Segments were the same when should have been different.");
        } catch (final IOException ioe) {
            Assert.fail("Unit test threw exception trying to read segments.", ioe);
        }
    }

    /**
     * Represent a Segment in the segmenter output.
     */
    private static class SegmentMean {
        public final double segmentMean;

        public SegmentMean(final double segmentMean) {
            this.segmentMean = segmentMean;
        }

        @Override
        public boolean equals(final Object other) {
            if (!(other instanceof SegmentMean))
                return false;
            final SegmentMean otherSegmentMean = (SegmentMean) other;
            return Math.abs(otherSegmentMean.segmentMean - segmentMean) < 0.0000000000001;
        }

        @Override
        public int hashCode() {
            return Double.hashCode(segmentMean);
        }
    }

    /**
     * Tsv reader for the Segmenter output.
     */
    private static class SegmentReader extends TableReader<SegmentMean> {
        public SegmentReader(final File file) throws IOException {
            super(file);
        }

        @Override
        protected SegmentMean createRecord(DataLine dataLine) {
            return new SegmentMean(dataLine.getDouble("Segment_Mean"));
        }
    }
}