package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Unit test for {@link BaselineCopyNumberCollection}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class BaselineCopyNumberCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/formats/collections");
    private static final File TEST_BASELINE_COPY_NUMBER_FILE = new File(TEST_SUB_DIR, "test_baseline_copy_number.tsv");
    private static final int[] EXPECTED_BASELINE_COPY_NUMBERS = new int[] {1, 2, 3, 4, 5, 0, 5, 4, 3, 2, 1};
    private static final String EXPECTED_SAMPLE_NAME = "TEST_SAMPLE_NAME";

    @Test
    public void testBaselineCopyNumberReading() {
        final BaselineCopyNumberCollection collection = new BaselineCopyNumberCollection(TEST_BASELINE_COPY_NUMBER_FILE);
        final List<IntegerCopyNumberState> expected = Arrays.stream(EXPECTED_BASELINE_COPY_NUMBERS)
                .mapToObj(IntegerCopyNumberState::new)
                .collect(Collectors.toList());
        Assert.assertEquals(collection.getRecords(), expected);
        Assert.assertEquals(collection.getMetadata().getSampleName(), EXPECTED_SAMPLE_NAME);
    }
}
