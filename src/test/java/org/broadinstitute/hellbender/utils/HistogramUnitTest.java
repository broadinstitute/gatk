package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.Histogram;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by gauthier on 11/1/16.
 */
public class HistogramUnitTest {
    private final Double EPSILON = 0.001;

    @Test
    public void testAdd() throws Exception {
        Histogram bimodalHist = new Histogram();
        for (int i = 0; i <= 100; i++) {
            bimodalHist.add(1 + i / 1000.0);
        }
        Assert.assertEquals(bimodalHist.get(1.0), new Integer(100), "");
        Assert.assertEquals(bimodalHist.get(1.1), new Integer(1), "");
    }

    @Test
    public void testAddingQuantizedValues() throws Exception {
        Histogram hist = new Histogram();
        for (int i = 0; i < 100; i++) {
            hist.add(1.2);
        }
        Assert.assertEquals(hist.get(1.2), new Integer(100));
        Assert.assertEquals(hist.median(), 1.2, EPSILON);
    }

    @Test
    public void testBulkAdd() throws Exception {
        Histogram bimodalHist = new Histogram();
        for (int i = 0; i <= 100; i++) {
            bimodalHist.add(1 + i / 1000.0, 2);
        }
        Assert.assertEquals(bimodalHist.get(1.0), new Integer(200), "");
        Assert.assertEquals(bimodalHist.get(1.1), new Integer(2), "");
    }

    @Test
    public void testMedianOfEvens() throws Exception {
        Histogram bimodalHist = new Histogram();
        for (int i = 0; i < 10; i++) {
            bimodalHist.add(10.0);
            bimodalHist.add(20.0);
        }

        Assert.assertEquals(bimodalHist.median(), 15.0, EPSILON, "");
    }

    @Test
    public void testMedianOfOdds() throws Exception {
        Histogram bimodalHist = new Histogram();
        for (int i = 0; i < 10; i++) {
            bimodalHist.add(10.0);
            bimodalHist.add(20.0);
        }
        bimodalHist.add(20.0);

        Assert.assertEquals(bimodalHist.median(), 20.0, EPSILON, "");
    }

    @Test
    public void testMedianOfEmptyHist() throws Exception {
        Histogram empty = new Histogram();
        Assert.assertNull(empty.median());
    }

    @Test
    public void testMedianOfSingleItem() throws Exception {
        Histogram singleItem = new Histogram();
        singleItem.add(20.0);
        Assert.assertEquals(singleItem.median(), 20.0, EPSILON, "");
    }
}

