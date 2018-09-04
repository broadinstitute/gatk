package org.broadinstitute.hellbender.tools.copynumber.utils.optimization;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Doubles;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link PersistenceOptimizer}.  Recovers local minima of simulated data,
 * e.g., data produced by adding noise to a sinusoid and a polynomial, along with other simple tests.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PersistenceOptimizerUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    @DataProvider(name = "dataPersistenceOptimizer")
    public Object[][] dataPersistenceOptimizer() {
        final int numPoints = 100;

        final Random rng = new Random(RANDOM_SEED);
        final double period = numPoints / 4.;
        final double[] dataSinusoidal = IntStream.range(0, numPoints)
                .mapToDouble(i -> Math.sin(i * 2 * Math.PI / period) + 0.1 * rng.nextGaussian())
                .toArray();
        final List<Integer> minimaIndicesExpectedSinusoidal = Arrays.asList(
                45, 94, 20, 68, 0, 82, 43, 22, 73, 31, 55, 79);                                         //from implementation at https://people.mpi-inf.mpg.de/~weinkauf/notes/persistence1d.html

        rng.setSeed(RANDOM_SEED);
        final double[] dataPolynomial = IntStream.range(0, numPoints)
                .mapToDouble(i -> (double) i / numPoints - 0.5)
                .map(x -> (x - 0.75) * (x - 0.5) * (x - 0.25) * (x - 0.1) * (x + 0.25) * (x + 0.5) + 0.0001 * rng.nextGaussian())
                .toArray();
        final List<Integer> minimaIndicesExpectedPolynomial = Arrays.asList(
                8, 68, 99, 73, 82, 89, 36, 34, 79, 45, 38, 49, 43, 84, 60, 56, 66, 62, 70, 92, 47, 87); //from implementation at https://people.mpi-inf.mpg.de/~weinkauf/notes/persistence1d.html

        final double[] dataIdentity = IntStream.range(0, numPoints)
                .mapToDouble(i -> i)
                .toArray();
        final List<Integer> minimaIndicesExpectedIdentity = Collections.singletonList(0);

        final double[] dataIdentityNegative = IntStream.range(0, numPoints)
                .mapToDouble(i -> -i)
                .toArray();
        final List<Integer> minimaIndicesExpectedIdentityNegative = Collections.singletonList(numPoints - 1);

        final double[] dataFlat = new double[numPoints];
        Arrays.fill(dataFlat, 1.);
        final List<Integer> minimaIndicesExpectedFlat = Collections.singletonList(0);               //leftmost point of a constant region is considered a local minimum

        final double[] dataFlatWithExtrema = new double[numPoints];
        Arrays.fill(dataFlatWithExtrema, 1.);
        dataFlatWithExtrema[10] = 0.;
        dataFlatWithExtrema[90] = 2.;
        final List<Integer> minimaIndicesExpectedFlatWithExtrema = Arrays.asList(10, 91, 0);        //leftmost point of a constant region is considered a local minimum

        final double[] dataAbsoluteValue = IntStream.range(0, numPoints)
                .mapToDouble(i -> Math.abs(i - numPoints / 2))
                .toArray();
        final List<Integer> minimaIndicesExpectedAbsoluteValue = Collections.singletonList(50);

        final double[] dataNegativeAbsoluteValue = IntStream.range(0, numPoints)
                .mapToDouble(i -> -Math.abs(i - numPoints / 2))
                .toArray();
        final List<Integer> minimaIndicesExpectedNegativeAbsoluteValue = Arrays.asList(0, 99);

        final double[] dataSingle = new double[]{1.};
        final List<Integer> minimaIndicesExpectedSingle = Collections.singletonList(0);

        return new Object[][]{
                {dataSinusoidal, minimaIndicesExpectedSinusoidal},
                {dataPolynomial, minimaIndicesExpectedPolynomial},
                {dataIdentity, minimaIndicesExpectedIdentity},
                {dataIdentityNegative, minimaIndicesExpectedIdentityNegative},
                {dataFlat, minimaIndicesExpectedFlat},
                {dataFlatWithExtrema, minimaIndicesExpectedFlatWithExtrema},
                {dataAbsoluteValue, minimaIndicesExpectedAbsoluteValue},
                {dataNegativeAbsoluteValue, minimaIndicesExpectedNegativeAbsoluteValue},
                {dataSingle, minimaIndicesExpectedSingle}
        };
    }

    @Test(dataProvider = "dataPersistenceOptimizer")
    public void testPersistenceOptimizer(final double[] data,
                                         final List<Integer> minimaIndicesExpected) {
        final PersistenceOptimizer optimizer = new PersistenceOptimizer(data);
        Assert.assertEquals(optimizer.getMinimaIndices(), minimaIndicesExpected);
        Assert.assertEquals(optimizer.getPersistences().get(0), Doubles.max(data) - Doubles.min(data));
        Assert.assertTrue(Ordering.natural().reverse().isOrdered(optimizer.getPersistences()));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyData() {
        new PersistenceOptimizer(new double[]{});
    }
}