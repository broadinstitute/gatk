package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;


public class SegmenterUnitTest extends BaseTest {

    @Test
    public void testHCC1143Reduced() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/segmenter/input/HCC1143_reduced.tsv");
        final File EXPECTED = new File("src/test/resources/segmenter/output/HCC1143_reduced_result.seg");
        final File output = createTempFile("recapseg.HCC1143", ".seg");
        String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testHCC1143Short() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/segmenter/input/HCC1143_short.tsv");
        final File EXPECTED = new File("src/test/resources/segmenter/output/HCC1143_short_result.seg");
        final File output = createTempFile("recapseg.HCC1143", ".seg");
        String sampleName = "HCC1143";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    @Test
    public void testSimple() throws IOException {
        final File INPUT_FILE = new File("src/test/resources/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/segmenter/output/Simple_result.seg");
        final File output = createTempFile("recapseg.HCC1143", ".seg");
        String sampleName = "Simple";
        RCBSSegmenter.writeSegmentFile(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), false);
        assertEqualSegments(output,EXPECTED);
    }

    /**
     * Compares the content of two segmenter output files.
     * @param actualOutput the actual segmenter output containing file.
     * @param expectedOutput the expected segmenter output containing file.
     * @throws NullPointerException if either {@code actualOutput} or {@code expectedOutput} is {@code null}.
     * @throws IOException if any was thrown when reading the input files.
     * @throws AssertionError if there are significant between both files.
     */
    public static void assertEqualSegments(final File actualOutput, final File expectedOutput) throws IOException {
        try (final SegmentReader actual = new SegmentReader(actualOutput);
             final SegmentReader expected = new SegmentReader(expectedOutput)) {
            final List<SegmentMean> actualSegmentMeans = actual.stream().collect(Collectors.toList());
            final List<SegmentMean> expectedSegmentMeans = expected.stream().collect(Collectors.toList());
            Assert.assertEquals(actualSegmentMeans, expectedSegmentMeans);
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