package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Unit test for {@link SimpleIntervalCollection}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class SimpleIntervalCollectionUnitTest extends GATKBaseTest {

    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/formats/collections");
    private static final SimpleIntervalCollection collection_0 = new SimpleIntervalCollection(
            new File(TEST_SUB_DIR, "interval_list_shard_0.tsv"));
    private static final SimpleIntervalCollection collection_1 = new SimpleIntervalCollection(
            new File(TEST_SUB_DIR, "interval_list_shard_1.tsv"));
    private static final SimpleIntervalCollection collection_2 = new SimpleIntervalCollection(
            new File(TEST_SUB_DIR, "interval_list_shard_2.tsv"));
    private static final SAMSequenceDictionary SAM_SEQUENCE_DICTIONARY = collection_0.getMetadata().getSequenceDictionary();

    @Test(dataProvider = "simpleIntervalCollectionSortedOrderTestData")
    public void testSimpleIntervalCollectionSortedOrder(final List<SimpleIntervalCollection> collections,
                                                        final int[] expectedOrder) {
        final int[] order = SimpleIntervalCollection.getSimpleIntervalCollectionSortedOrder(
                collections, SAM_SEQUENCE_DICTIONARY, true);
        ArrayAsserts.assertArrayEquals(expectedOrder, order);
    }

    @Test(dataProvider = "simpleIntervalCollectionSortedOrderBadTestData", expectedExceptions = IllegalArgumentException.class)
    public void testSimpleIntervalCollectionSortedOrderBadInput(final List<SimpleIntervalCollection> collections) {
        SimpleIntervalCollection.getSimpleIntervalCollectionSortedOrder(collections, SAM_SEQUENCE_DICTIONARY, true);
    }

    @DataProvider(name = "simpleIntervalCollectionSortedOrderTestData")
    public Object[][] getSimpleIntervalCollectionSortedOrderTestData() {
        return new Object[][] {
                {Arrays.asList(collection_0, collection_1, collection_2), new int[] {0, 1, 2}},
                {Arrays.asList(collection_1, collection_0, collection_2), new int[] {1, 0, 2}},
                {Arrays.asList(collection_2, collection_0, collection_1), new int[] {1, 2, 0}}
        };
    }

    @DataProvider(name = "simpleIntervalCollectionSortedOrderBadTestData")
    public Object[][] getSimpleIntervalCollectionSortedOrderBadTestData() {
        return new Object[][] { /* repeated collections */
                {Arrays.asList(collection_0, collection_0, collection_2)},
                {Arrays.asList(collection_2, collection_2, collection_2)},
                {Arrays.asList(collection_2, collection_1, collection_1)}
        };
    }
}
