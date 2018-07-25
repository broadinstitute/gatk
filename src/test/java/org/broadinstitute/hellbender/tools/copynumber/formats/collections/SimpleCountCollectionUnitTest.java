package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


public final class SimpleCountCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/formats/collections");
    private static final File INTEGER_COUNTS_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts.tsv");
    private static final File INTEGER_COUNTS_HDF5_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts.hdf5");
    private static final File INTEGER_COUNTS_MISSING_HEADER_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts-missing-header.tsv");
    private static final File DOUBLE_COUNTS_FILE = new File(TEST_SUB_DIR, "simple-count-collection-double-counts.tsv");

    private static final SampleLocatableMetadata METADATA_EXPECTED = new SimpleSampleLocatableMetadata(
            "test-sample",
            new SAMSequenceDictionary(Collections.singletonList(
                    new SAMSequenceRecord("20", 200000))));

    private static final List<SimpleInterval> INTERVALS_EXPECTED = Arrays.asList(
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
    private static final RealMatrix READ_COUNTS_EXPECTED = new Array2DRowRealMatrix(
            new double[][]{{0, 0, 0, 0, 0, 0, 94, 210, 22, 21, 24, 84, 247, 181, 27, 72}});

    @Test
    public void testReadIntegerCounts() {
        final SimpleCountCollection scc = SimpleCountCollection.read(INTEGER_COUNTS_FILE);
        final SampleLocatableMetadata metadata = scc.getMetadata();
        final List<SimpleInterval> intervals = scc.getIntervals();
        final RealMatrix readCounts = new Array2DRowRealMatrix(new double[][]{scc.getRecords().stream().mapToDouble(SimpleCount::getCount).toArray()});

        Assert.assertEquals(metadata, METADATA_EXPECTED);
        Assert.assertEquals(intervals, INTERVALS_EXPECTED);
        Assert.assertEquals(readCounts, READ_COUNTS_EXPECTED);
    }

    @Test
    public void testReadIntegerCountsHDF5() {
        final SimpleCountCollection scc = SimpleCountCollection.read(INTEGER_COUNTS_HDF5_FILE);
        final SampleLocatableMetadata metadata = scc.getMetadata();
        final List<SimpleInterval> intervals = scc.getIntervals();
        final RealMatrix readCounts = new Array2DRowRealMatrix(new double[][]{scc.getRecords().stream().mapToDouble(SimpleCount::getCount).toArray()});

        Assert.assertEquals(metadata, METADATA_EXPECTED);
        Assert.assertEquals(intervals, INTERVALS_EXPECTED);
        Assert.assertEquals(readCounts, READ_COUNTS_EXPECTED);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReadIntegerCountsMissingHeader() {
        SimpleCountCollection.read(INTEGER_COUNTS_MISSING_HEADER_FILE);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadDoubleCounts() {
        SimpleCountCollection.read(DOUBLE_COUNTS_FILE);
    }
}