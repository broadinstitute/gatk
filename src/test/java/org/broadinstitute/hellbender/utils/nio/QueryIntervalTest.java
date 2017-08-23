package org.broadinstitute.hellbender.utils.nio;

import htsjdk.samtools.QueryInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Test for QueryInterval
 */
public class QueryIntervalTest {

  @Test
  public void testAbuts() throws Exception {
    QueryInterval int0 = new QueryInterval(0, 1519, 1520);
    // overlaps
    QueryInterval int1 = new QueryInterval(0, 1520, 1522);
    // abuts
    QueryInterval int2 = new QueryInterval(0, 1521, 1522);
    // separated
    QueryInterval int3 = new QueryInterval(0, 1522, 1522);

    Assert.assertFalse(int0.abuts(int1));
    Assert.assertTrue(int0.abuts(int2));
    Assert.assertFalse(int0.abuts(int3));
  }

  @Test
  public void testOverlaps() throws Exception {
  }

  @Test
  public void testOptimizeIntervals() throws Exception {
    QueryInterval[] intervals = new QueryInterval[] {
        new QueryInterval(0, 1519, 1520),
        new QueryInterval(0, 1521, 1525)
    };
    QueryInterval[] expected = new QueryInterval[] {
        new QueryInterval(0, 1519, 1525),
    };
    final QueryInterval[] got = QueryInterval.optimizeIntervals(intervals);
    Assert.assertEquals(expected, got);
  }

}