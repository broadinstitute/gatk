package org.broadinstitute.hellbender.utils;


// the imports for unit testing.


import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public final class MedianUnitTest extends BaseTest {

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private class MedianTestProvider extends TestDataProvider {
        final List<Integer> values = new ArrayList<>();
        final int cap;
        final Integer expected;

        public MedianTestProvider(int expected, int cap, Integer... values) {
            super(MedianTestProvider.class);
            this.expected = expected;
            this.cap = cap;
            this.values.addAll(Arrays.asList(values));
            this.name = String.format("values=%s expected=%d cap=%d", this.values, expected, cap);
        }
    }

    @DataProvider(name = "MedianTestProvider")
    public Object[][] makeMedianTestProvider() {
        new MedianTestProvider(1, 1000, 0, 1, 2);
        new MedianTestProvider(1, 1000, 1, 0, 1, 2);
        new MedianTestProvider(1, 1000, 0, 1, 2, 3);
        new MedianTestProvider(2, 1000, 0, 1, 2, 3, 4);
        new MedianTestProvider(2, 1000, 4, 1, 2, 3, 0);
        new MedianTestProvider(1, 1000, 1);
        new MedianTestProvider(2, 1000, 2);
        new MedianTestProvider(1, 1000, 1, 2);

        new MedianTestProvider(1, 3, 1);
        new MedianTestProvider(1, 3, 1, 2);
        new MedianTestProvider(2, 3, 1, 2, 3);
        new MedianTestProvider(2, 3, 1, 2, 3, 4);
        new MedianTestProvider(2, 3, 1, 2, 3, 4, 5);

        new MedianTestProvider(1, 3, 1);
        new MedianTestProvider(1, 3, 1, 2);
        new MedianTestProvider(2, 3, 3, 2, 1);
        new MedianTestProvider(3, 3, 4, 3, 2, 1);
        new MedianTestProvider(4, 3, 5, 4, 3, 2, 1);

        return MedianTestProvider.getTests(MedianTestProvider.class);
    }

    @Test(dataProvider = "MedianTestProvider")
    public void testBasicLikelihoods(MedianTestProvider cfg) {
        final Median<Integer> median = new Median<>(cfg.cap);

        int nAdded = 0;
        for ( final int value : cfg.values )
            if ( median.add(value) )
                nAdded++;

        Assert.assertEquals(nAdded, median.size());

        Assert.assertEquals(cfg.values.isEmpty(), median.isEmpty());
        Assert.assertEquals(cfg.values.size() >= cfg.cap, median.isFull());
        Assert.assertEquals(median.getMedian(), cfg.expected, cfg.toString());
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testEmptyMedian() {
        final Median<Integer> median = new Median<>();
        Assert.assertTrue(median.isEmpty());
        final Integer d = 100;
        Assert.assertEquals(median.getMedian(d), d);
        median.getMedian();
    }

}