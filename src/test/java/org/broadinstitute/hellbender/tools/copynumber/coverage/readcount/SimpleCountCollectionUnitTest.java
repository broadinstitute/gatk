package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public final class SimpleCountCollectionUnitTest extends GATKBaseTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "copynumber/coverage/readcount/";
    private static final File INTEGER_COUNTS_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts.tsv");
    private static final File DOUBLE_COUNTS_FILE = new File(TEST_SUB_DIR, "simple-count-collection-double-counts.tsv");

    private static final String SAMPLE_NAME_EXPECTED = "test";
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
        final String sampleName = scc.getSampleName();
        final List<SimpleInterval> intervals = scc.getIntervals();
        final RealMatrix readCounts = new Array2DRowRealMatrix(new double[][]{scc.getRecords().stream().mapToDouble(SimpleCount::getCount).toArray()});

        Assert.assertEquals(sampleName, SAMPLE_NAME_EXPECTED);
        Assert.assertEquals(intervals, INTERVALS_EXPECTED);
        Assert.assertEquals(readCounts, READ_COUNTS_EXPECTED);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadDoubleCounts() {
        SimpleCountCollection.read(DOUBLE_COUNTS_FILE);
    }
}