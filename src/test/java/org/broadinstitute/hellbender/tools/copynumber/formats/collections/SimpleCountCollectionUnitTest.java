package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public final class SimpleCountCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/formats/collections");

    private static final File INTEGER_COUNTS_TSV_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts.tsv");
    private static final File INTEGER_COUNTS_HDF5_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts.hdf5");
    private static final File INTEGER_COUNTS_TSV_GZ_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts.tsv.gz");

    private static final String GCS_COLLECTIONS_TEST_SUBDIRECTORY = "org/broadinstitute/hellbender/tools/copynumber/formats/collections/";
    private static final String GCS_INTEGER_COUNTS_TSV_SUB_PATH = GCS_COLLECTIONS_TEST_SUBDIRECTORY +
            "simple-count-collection-integer-counts.counts.tsv";
    private static final String GCS_INTEGER_COUNTS_TSV_GZ_SUB_PATH = GCS_COLLECTIONS_TEST_SUBDIRECTORY +
            "simple-count-collection-integer-counts.counts.tsv.gz";

    private static final File INTEGER_COUNTS_MISSING_HEADER_TSV_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts-missing-header.tsv");
    private static final File DOUBLE_COUNTS_TSV_FILE = new File(TEST_SUB_DIR, "simple-count-collection-double-counts.tsv");

    private static final SampleLocatableMetadata METADATA_EXPECTED = new SimpleSampleLocatableMetadata(
            "test-sample",
            new SAMSequenceDictionary(Collections.singletonList(
                    new SAMSequenceRecord("20", 200000))));
    private static final List<SimpleInterval> INTERVALS_EXPECTED = ImmutableList.of(
            new SimpleInterval("20", 1,10000),
            new SimpleInterval("20", 10001,20000),
            new SimpleInterval("20", 20001, 30000),
            new SimpleInterval("20", 30001, 40000),
            new SimpleInterval("20", 40001, 50000),
            new SimpleInterval("20", 50001, 60000),
            new SimpleInterval("20", 60001, 70000),
            new SimpleInterval("20", 70001, 80000),
            new SimpleInterval("20", 80001, 90000),
            new SimpleInterval("20", 90001, 100000),
            new SimpleInterval("20", 100001, 110000),
            new SimpleInterval("20", 110001, 120000),
            new SimpleInterval("20", 120001, 130000),
            new SimpleInterval("20", 130001, 140000),
            new SimpleInterval("20", 140001, 150000),
            new SimpleInterval("20", 150001, 160000));
    private static final List<Integer> COUNTS_EXPECTED =
            ImmutableList.of(0, 0, 0, 0, 0, 0, 94, 210, 22, 21, 24, 84, 247, 181, 27, 72);
    private static final List<SimpleCount> SIMPLE_COUNTS_EXPECTED = IntStream.range(0, INTERVALS_EXPECTED.size()).boxed()
            .map(i -> new SimpleCount(INTERVALS_EXPECTED.get(i), COUNTS_EXPECTED.get(i)))
            .collect(Collectors.toList());


    private static void assertCountsExpected(final SimpleCountCollection counts) {
        Assert.assertEquals(counts.getMetadata(), METADATA_EXPECTED);
        Assert.assertEquals(counts.getIntervals(), INTERVALS_EXPECTED);
        Assert.assertEquals(counts.getRecords(), SIMPLE_COUNTS_EXPECTED);
    }

    private static void assertCountsSubsetExpected(final SimpleCountCollection countsSubset,
                                                   final List<Integer> indices) {
        Assert.assertEquals(countsSubset.getMetadata(), METADATA_EXPECTED);
        Assert.assertEquals(countsSubset.getIntervals(), indices.stream().map(INTERVALS_EXPECTED::get).collect(Collectors.toList()));
        Assert.assertEquals(countsSubset.getRecords(), indices.stream().map(SIMPLE_COUNTS_EXPECTED::get).collect(Collectors.toList()));
    }

    @DataProvider(name = "simpleCountCollectionReadTestData")
    public Object[] getSimpleCountCollectionReadTestData() {
        return new Object[] {
                INTEGER_COUNTS_TSV_FILE,
                INTEGER_COUNTS_HDF5_FILE,
                INTEGER_COUNTS_TSV_GZ_FILE
        };
    }

    @Test(dataProvider = "simpleCountCollectionReadTestData")
    public void testRead(final File file) {
        final SimpleCountCollection counts = SimpleCountCollection.read(file);
        assertCountsExpected(counts);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReadMissingHeader() {
        SimpleCountCollection.read(INTEGER_COUNTS_MISSING_HEADER_TSV_FILE);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadDoubleCounts() {
        SimpleCountCollection.read(DOUBLE_COUNTS_TSV_FILE);
    }

    @DataProvider(name = "simpleCountCollectionReadFromGCSTestData")
    public Object[] getSimpleCountCollectionReadFromGCSTestData() {
        return new Object[] {
                GCS_INTEGER_COUNTS_TSV_SUB_PATH,
                GCS_INTEGER_COUNTS_TSV_GZ_SUB_PATH
        };
    }

    @Test(dataProvider = "simpleCountCollectionReadFromGCSTestData", groups = {"bucket"})
    public void testReadFromGCS(final String subPath) {
        final SimpleCountCollection counts = SimpleCountCollection.readFromGCS(GCloudTestUtils.getTestInputPath() + subPath);
        assertCountsExpected(counts);
    }

    @Test(dataProvider = "simpleCountCollectionReadTestData")
    public void testReadAndSubset(final File file) {
        final SimpleCountCollection countsSubsetNull = SimpleCountCollection.readAndSubset(file,
                null);
        assertCountsExpected(countsSubsetNull);

        final SimpleCountCollection countsSubsetEmpty = SimpleCountCollection.readAndSubset(file,
                Collections.emptySet());
        assertCountsExpected(countsSubsetEmpty);

        final SimpleCountCollection countsSubset1To9 = SimpleCountCollection.readAndSubset(file,
                new HashSet<>(INTERVALS_EXPECTED.subList(1, 10)));
        assertCountsSubsetExpected(countsSubset1To9, Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9));

        final SimpleCountCollection countsSubsetFirstThreeIntervalsOutOfOrder = SimpleCountCollection.readAndSubset(file,
                ImmutableSet.of(
                        new SimpleInterval("20", 10001,20000),
                        new SimpleInterval("20", 1,10000),
                        new SimpleInterval("20", 20001,30000)));
        assertCountsSubsetExpected(countsSubsetFirstThreeIntervalsOutOfOrder, Arrays.asList(0, 1, 2));  //first three intervals

        final SimpleCountCollection countsSubsetThreeIntervals = SimpleCountCollection.readAndSubset(file,
                ImmutableSet.of(
                        new SimpleInterval("20", 10001,20000),
                        new SimpleInterval("20", 70001,80000),
                        new SimpleInterval("20", 140001,150000)));
        assertCountsSubsetExpected(countsSubsetThreeIntervals, Arrays.asList(1, 7, 14));                //three non-adjacent intervals

        final SimpleCountCollection countsSubsetFirstAndAHalfIntervals = SimpleCountCollection.readAndSubset(file,
                ImmutableSet.of(
                        new SimpleInterval("20", 1, 10000),
                        new SimpleInterval("20", 10001, 15000)));
        assertCountsSubsetExpected(countsSubsetFirstAndAHalfIntervals, Collections.singletonList(0));   //first interval

        final SimpleCountCollection countsSubsetEntireRangeSingleInterval = SimpleCountCollection.readAndSubset(file,
                new HashSet<>(Collections.singletonList(new SimpleInterval("20", 1, 160000))));
        assertCountsSubsetExpected(countsSubsetEntireRangeSingleInterval, Collections.emptyList());      //empty list
    }

    @Test(dataProvider = "simpleCountCollectionReadFromGCSTestData", groups = {"bucket"})
    public void testReadOverlappingSubsetFromGCS(final String subPath) {
        final String path = GCloudTestUtils.getTestInputPath() + subPath;

        final SimpleCountCollection countsSubsetNull = SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                null);
        assertCountsExpected(countsSubsetNull);

        final SimpleCountCollection countsSubsetEmpty = SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Collections.emptyList());
        assertCountsExpected(countsSubsetEmpty);

        final SimpleCountCollection countsSubsetFirstThreeIntervalsExact = SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Arrays.asList(
                        new SimpleInterval("20", 1, 10000),
                        new SimpleInterval("20", 10001, 20000),
                        new SimpleInterval("20", 20001, 30000)));
        assertCountsSubsetExpected(countsSubsetFirstThreeIntervalsExact, Arrays.asList(0, 1, 2));           //first three intervals

        final SimpleCountCollection countsSubsetFirstThreeIntervalsExactRange = SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Collections.singletonList(new SimpleInterval("20", 1, 30000)));
        assertCountsSubsetExpected(countsSubsetFirstThreeIntervalsExactRange, Arrays.asList(0, 1, 2));      //first three intervals

        final SimpleCountCollection countsSubsetFirstThreeIntervalsOverlapping = SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Collections.singletonList(new SimpleInterval("20", 5000, 20001)));
        assertCountsSubsetExpected(countsSubsetFirstThreeIntervalsOverlapping, Arrays.asList(0, 1, 2));     //first three intervals

        final SimpleCountCollection countsSubsetThreeIntervalsOverlapping = SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Arrays.asList(
                        new SimpleInterval("20", 1, 15000),
                        new SimpleInterval("20", 145000, 150000)));
        assertCountsSubsetExpected(countsSubsetThreeIntervalsOverlapping, Arrays.asList(0, 1, 14));         //three intervals

        final SimpleCountCollection countsSubsetEntireRangeTwoIntervals = SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Arrays.asList(
                        new SimpleInterval("20", 1, 50000),
                        new SimpleInterval("20", 50001, 160000)));
        assertCountsExpected(countsSubsetEntireRangeTwoIntervals);                                          //entire list

        final SimpleCountCollection countsSubsetEntireRangeOneInterval = SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Collections.singletonList(new SimpleInterval("20", 1, 160000)));
        assertCountsExpected(countsSubsetEntireRangeOneInterval);                                           //entire list
    }

    @Test(dataProvider = "simpleCountCollectionReadFromGCSTestData", groups = {"bucket"}, expectedExceptions = IllegalArgumentException.class)
    public void testReadOverlappingSubsetFromGCSSortedIntervals(final String subPath) {
        final String path = GCloudTestUtils.getTestInputPath() + subPath;

        SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Arrays.asList(
                        new SimpleInterval("20", 10001,20000),
                        new SimpleInterval("20", 1,10000),
                        new SimpleInterval("20", 20001,30000)));
    }

    @Test(dataProvider = "simpleCountCollectionReadFromGCSTestData", groups = {"bucket"}, expectedExceptions = IllegalArgumentException.class)
    public void testReadOverlappingSubsetFromGCSNonOverlappingIntervals(final String subPath) {
        final String path = GCloudTestUtils.getTestInputPath() + subPath;

        SimpleCountCollection.readOverlappingSubsetFromGCS(path,
                Arrays.asList(
                        new SimpleInterval("20", 1,20000),
                        new SimpleInterval("20", 10000,20000)));
    }
}
