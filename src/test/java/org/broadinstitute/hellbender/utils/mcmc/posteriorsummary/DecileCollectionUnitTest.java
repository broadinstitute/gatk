package org.broadinstitute.hellbender.utils.mcmc.posteriorsummary;

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
    private static final List<Double> testDeciles = IntStream.range(0, 11).boxed().map(Integer::doubleValue).collect(Collectors.toList());

    @Test
    public void testConstructorSamplesMode() {
        final DecileCollection result = new DecileCollection(testSamples, DecileCollection.ConstructionMode.SAMPLES);
        Assert.assertEquals(result.get(Decile.DECILE_0), 1.);
        Assert.assertEquals(result.get(Decile.DECILE_10), 10.);
        Assert.assertEquals(result.get(Decile.DECILE_20), 20.);
        Assert.assertEquals(result.get(Decile.DECILE_30), 30.);
        Assert.assertEquals(result.get(Decile.DECILE_40), 40.);
        Assert.assertEquals(result.get(Decile.DECILE_50), 50.);
        Assert.assertEquals(result.get(Decile.DECILE_60), 60.);
        Assert.assertEquals(result.get(Decile.DECILE_70), 70.);
        Assert.assertEquals(result.get(Decile.DECILE_80), 80.);
        Assert.assertEquals(result.get(Decile.DECILE_90), 90.);
        Assert.assertEquals(result.get(Decile.DECILE_100), 99.);
    }

    @Test
    public void testConstructorDecilesMode() {
        final DecileCollection result = new DecileCollection(testDeciles, DecileCollection.ConstructionMode.DECILES);
        Assert.assertEquals(result.get(Decile.DECILE_0), 0.);
        Assert.assertEquals(result.get(Decile.DECILE_10), 1.);
        Assert.assertEquals(result.get(Decile.DECILE_20), 2.);
        Assert.assertEquals(result.get(Decile.DECILE_30), 3.);
        Assert.assertEquals(result.get(Decile.DECILE_40), 4.);
        Assert.assertEquals(result.get(Decile.DECILE_50), 5.);
        Assert.assertEquals(result.get(Decile.DECILE_60), 6.);
        Assert.assertEquals(result.get(Decile.DECILE_70), 7.);
        Assert.assertEquals(result.get(Decile.DECILE_80), 8.);
        Assert.assertEquals(result.get(Decile.DECILE_90), 9.);
        Assert.assertEquals(result.get(Decile.DECILE_100), 10.);
    }

    @Test
    public void testConstructorNaN() {
        final DecileCollection result = new DecileCollection(Collections.singletonList(Double.NaN), DecileCollection.ConstructionMode.SAMPLES);
        Assert.assertEquals(result.getAll(), Collections.nCopies(DecileCollection.NUM_DECILES, Double.NaN));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmpty() {
        new DecileCollection(Collections.emptyList(), DecileCollection.ConstructionMode.SAMPLES);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDecilesModeIncorrectNumberOfValues() {
        new DecileCollection(Arrays.asList(0., 1., 2.), DecileCollection.ConstructionMode.DECILES);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDecilesModeUnsortedValues() {
        new DecileCollection(Arrays.asList(0., 1., 2., 3., 4., 6., 5., 7., 8., 9., 10.), DecileCollection.ConstructionMode.DECILES);
    }

    @Test
    public void testGetAll() {
        final DecileCollection result = new DecileCollection(testSamples, DecileCollection.ConstructionMode.SAMPLES);
        Assert.assertEquals(result.getAll(), Arrays.asList(1., 10., 20., 30., 40., 50., 60., 70., 80., 90., 99.));
    }

    @Test
    public void testGetInner() {
        final DecileCollection result = new DecileCollection(testSamples, DecileCollection.ConstructionMode.SAMPLES);
        Assert.assertEquals(result.getInner(), Arrays.asList(10., 20., 30., 40., 50., 60., 70., 80., 90.));
    }
}