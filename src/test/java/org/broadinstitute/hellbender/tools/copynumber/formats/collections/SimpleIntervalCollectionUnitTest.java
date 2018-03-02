package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
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

    @Test(dataProvider = "simpleIntervalCollectionSortedOrderTestData")
    public void testSimpleIntervalCollectionSortedOrder(final List<SimpleIntervalCollection> collections,
                                                        final List<Integer> expectedOrder) {
        final List<Integer> order = AbstractLocatableCollection.getShardedCollectionSortOrder(collections);
        Assert.assertEquals(expectedOrder, order);
    }

    @Test(dataProvider = "simpleIntervalCollectionSortedOrderBadTestData", expectedExceptions = IllegalArgumentException.class)
    public void testSimpleIntervalCollectionSortedOrderBadInput(final List<SimpleIntervalCollection> collections) {
        SimpleIntervalCollection.getShardedCollectionSortOrder(collections);
    }

    @DataProvider(name = "simpleIntervalCollectionSortedOrderTestData")
    public Object[][] getSimpleIntervalCollectionSortedOrderTestData() {
        return new Object[][] {
                {Arrays.asList(collection_0, collection_1, collection_2), Arrays.asList(0, 1, 2)},
                {Arrays.asList(collection_1, collection_0, collection_2), Arrays.asList(1, 0, 2)},
                {Arrays.asList(collection_2, collection_0, collection_1), Arrays.asList(1, 2, 0)}
        };
    }

    @DataProvider(name = "simpleIntervalCollectionSortedOrderBadTestData")
    public Object[][] getSimpleIntervalCollectionSortedOrderBadTestData() {
        final List<SAMSequenceRecord> extendedSAMSequenceRecords = new ArrayList<>(
                collection_0.getMetadata().getSequenceDictionary().getSequences());
        extendedSAMSequenceRecords.add(new SAMSequenceRecord("a_very_special_contig", 10));
        final SAMSequenceDictionary extendedSAMSequenceDictionary = new SAMSequenceDictionary(extendedSAMSequenceRecords);
        final SimpleIntervalCollection collection_0_bad = new SimpleIntervalCollection(
                new SimpleLocatableMetadata(extendedSAMSequenceDictionary),
                collection_0.getRecords());
        return new Object[][] { /* repeated collections */
                {Arrays.asList(collection_0, collection_0, collection_2)},
                {Arrays.asList(collection_2, collection_2, collection_2)},
                {Arrays.asList(collection_2, collection_1, collection_1)},
                {Arrays.asList(collection_0_bad, collection_1, collection_2)} /* collect_0_bad has different metadata */
        };
    }
}
