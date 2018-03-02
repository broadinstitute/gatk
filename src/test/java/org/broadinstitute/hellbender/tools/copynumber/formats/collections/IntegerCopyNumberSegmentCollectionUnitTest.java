package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Unit test for {@link IntegerCopyNumberSegmentCollection}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberSegmentCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/formats/collections");
    public static final File TEST_INTEGER_COPY_NUMBER_SEGMENTS_FILE =
            new File(TEST_SUB_DIR, "test_copy_number_segments.tsv");
    public static final String EXPECTED_SAMPLE_NAME = "TEST_SAMPLE_NAME";

    /* the following segments are expected to be in {@link #TEST_INTEGER_COPY_NUMBER_SEGMENTS_FILE} */
    public static final List<IntegerCopyNumberSegment> TEST_INTEGER_COPY_NUMBER_SEGMENTS = Arrays.asList(
            new IntegerCopyNumberSegment(
                    new SimpleInterval("1", 30365, 30503),
                    new IntegerCopyNumberState(2),
                    new IntegerCopyNumberState(3),
                    1, 115.84, 115.84, 115.84, 115.84),
            new IntegerCopyNumberSegment(
                    new SimpleInterval("2", 41605, 224922),
                    new IntegerCopyNumberState(3),
                    new IntegerCopyNumberState(2),
                    4, 49.14, 51.24, 101.62, 51.50),
            new IntegerCopyNumberSegment(
                    new SimpleInterval("3", 1189690, 1262499),
                    new IntegerCopyNumberState(4),
                    new IntegerCopyNumberState(2),
                    2, 53.83, 28.88, 28.88, 104.34),
            new IntegerCopyNumberSegment(
                    new SimpleInterval("X", 2986139, 7177815),
                    new IntegerCopyNumberState(5),
                    new IntegerCopyNumberState(1),
                    41, 192.92, 192.92, 212.98, 250.44),
            new IntegerCopyNumberSegment(
                    new SimpleInterval("Y", 150854, 249631),
                    new IntegerCopyNumberState(6),
                    new IntegerCopyNumberState(1),
                    15, 195.77, 195.77, 266.47, 311.75));

    @Test
    public void testIntegerCopyNumberSegmentCollectionReading() {
        final IntegerCopyNumberSegmentCollection collection = new IntegerCopyNumberSegmentCollection(
                TEST_INTEGER_COPY_NUMBER_SEGMENTS_FILE);
        Assert.assertEquals(collection.getRecords(), TEST_INTEGER_COPY_NUMBER_SEGMENTS);
        Assert.assertEquals(collection.getMetadata().getSampleName(), EXPECTED_SAMPLE_NAME);
    }

    @Test
    public void testWriteRead() {
        final SampleLocatableMetadata metadata = new IntegerCopyNumberSegmentCollection(
                TEST_INTEGER_COPY_NUMBER_SEGMENTS_FILE).getMetadata();
        final IntegerCopyNumberSegmentCollection collection = new IntegerCopyNumberSegmentCollection(
                metadata, TEST_INTEGER_COPY_NUMBER_SEGMENTS);
        final File tempCollectionFile = createTempFile("integer-copy-number-segments-collection", ".tsv");
        collection.write(tempCollectionFile);
        final IntegerCopyNumberSegmentCollection readCollection = new IntegerCopyNumberSegmentCollection(
                tempCollectionFile);
        Assert.assertEquals(readCollection.getMetadata().getSampleName(), EXPECTED_SAMPLE_NAME);
        Assert.assertEquals(readCollection.getMetadata(), metadata);
        Assert.assertEquals(readCollection.getRecords(), TEST_INTEGER_COPY_NUMBER_SEGMENTS);
    }
}
