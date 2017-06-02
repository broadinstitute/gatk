package org.broadinstitute.hellbender.utils.mcmc;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tests for {@link DecileCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class DecileCollectionUnitTest {

    private static final List<Double> testSamples = IntStream.range(1, 100).boxed().map(Integer::doubleValue).collect(Collectors.toList());

    @Test
    public void testConstructor() {
        final DecileCollection result = new DecileCollection(testSamples);
        Assert.assertEquals(result.get(Decile.DECILE_10), 10.);
        Assert.assertEquals(result.get(Decile.DECILE_20), 20.);
        Assert.assertEquals(result.get(Decile.DECILE_30), 30.);
        Assert.assertEquals(result.get(Decile.DECILE_40), 40.);
        Assert.assertEquals(result.get(Decile.DECILE_50), 50.);
        Assert.assertEquals(result.get(Decile.DECILE_60), 60.);
        Assert.assertEquals(result.get(Decile.DECILE_70), 70.);
        Assert.assertEquals(result.get(Decile.DECILE_80), 80.);
        Assert.assertEquals(result.get(Decile.DECILE_90), 90.);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmpty() {
        new DecileCollection(Collections.emptyList());
    }

    @Test
    public void testGetAll() {
        final DecileCollection result = new DecileCollection(testSamples);
        Assert.assertEquals(result.getAll(), Arrays.asList(10., 20., 30., 40., 50., 60., 70., 80., 90.));
    }
}