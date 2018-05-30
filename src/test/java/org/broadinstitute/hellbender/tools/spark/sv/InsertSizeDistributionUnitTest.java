package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.IntToDoubleFunction;

/**
 * Unit tests for {@link InsertSizeDistribution}.
 */
public class InsertSizeDistributionUnitTest extends BaseTest {

    @Test(dataProvider = "testData")
    public void testProbability(final String description, final int x, final double expected, final double logExpected) {
        final InsertSizeDistribution isd = new InsertSizeDistribution(description);
            Assert.assertEquals(isd.density(x), expected, 0.00001);
        Assert.assertEquals(isd.logDensity(x), logExpected, 0.01);
    }

    @DataProvider(name = "testData")
    public Object[][] testData() {
        final List<Object[]> result = new ArrayList<>();
        // Expected values calculated using R.
        // dnorm(231, 300, 150) == 0.002392602
        final Random random = new Random(131);
        final JDKRandomGenerator randomGenerator = new JDKRandomGenerator();
        randomGenerator.setSeed(random.nextInt());
        final double[] means = {100, 200, 300, 10.5, 310.12, 10313.0};
        final double[] cvs = {0.01, 0.1, 0.5, 1, 2};
        final int[] fixedSizes = {1, 11, 113, 143, 243, 321, 494, 539, 10190, 301298712};
        final double[] sizeSigmas = {0, -1, 1, 3.5, -2, 2, -6.9, 6.9};
        final PoissonDistribution spacesDistr = new PoissonDistribution(randomGenerator, 0.1, 0.0001, 100);
        for (final InsertSizeDistribution.Type type : InsertSizeDistribution.SUPPORTED_TYPES) {
            final List<String> distrNames = type.getNames();
            for (final double mean : means) {
                for (final double cv : cvs) {
                    final double stddev = mean * cv;
                    // We use alternative code to compose the densities using the formulas found in wikipedia articles.
                    // the actual implementation in main relies on apache common math, so is difficult to fall into the
                    // same error mode thus masking bugs.
                    final IntToDoubleFunction expectedDensity;
                    final IntToDoubleFunction logExpecetedDensity;
                    if (type.getClass() == InsertSizeDistribution.NormalType.class) {
                        expectedDensity = (x) -> Math.exp(-.5 * Math.pow((((double) x) - mean) / stddev, 2)) * (1.0 / (stddev * Math.sqrt(2 * Math.PI)));
                        logExpecetedDensity = (x) -> -.5 * Math.pow((((double) x) - mean) / stddev, 2) - Math.log(stddev * Math.sqrt(2 * Math.PI));
                    } else if (type.getClass() == InsertSizeDistribution.LogNormalType.class) {
                        final double var = stddev * stddev;
                        final double logMean = Math.log(mean) - Math.log(Math.sqrt(1 + (var / (mean * mean))));
                        final double logStddev = Math.sqrt(Math.log(1 + var / (mean * mean)));
                        expectedDensity = (x) -> Math.exp(-.5 * Math.pow((Math.log(x) - logMean) / logStddev, 2)) / (x * logStddev * Math.sqrt(2 * Math.PI));
                        logExpecetedDensity = (x) -> -.5 * Math.pow((Math.log(x) - logMean) / logStddev, 2) - Math.log(x * logStddev * Math.sqrt(2 * Math.PI));
                    } else {
                        throw new IllegalStateException("test do not support one of the type supported by InsertSizeDistribution: " + type.getNames().get(0));
                    }
                    // We add fixed length cases
                    for (final int fixedSize : fixedSizes) {
                        final String distrName = distrNames.get(random.nextInt(distrNames.size()));
                        result.add(new Object[]{
                                composeDescriptionString(distrName, mean, stddev, spacesDistr),
                                fixedSize, expectedDensity.applyAsDouble(fixedSize),
                        logExpecetedDensity.applyAsDouble(fixedSize)});
                    }
                    // We add relative length cases (expressed in sigmas)
                    for (final double sizeSigma : sizeSigmas) {
                        final int x = (int) Math.round(sizeSigma * stddev + mean);
                        if (x <= 0) {
                            continue; // we skip non-sense sizes.
                        }
                        final String distrName = distrNames.get(random.nextInt(distrNames.size()));
                        result.add(new Object[]{
                                composeDescriptionString(distrName, mean, stddev, spacesDistr),
                                x, expectedDensity.applyAsDouble(x),
                        logExpecetedDensity.applyAsDouble(x)});
                    }
                }
            }
        }
        // A couple of hard-wired cases to attest that the code above is generating genuine cases
        // rather than have the same error mode as the implementation in main.
        result.add(new Object[] { "N(300,150)", 231, 0.002392602, Math.log(0.002392602)});
        result.add(new Object[] { "lnN(100,1)", 103, 0.004833924, Math.log(0.004833924)});
        return result.toArray(new Object[result.size()][]);
    }

    /**
     * Composes a isd descriptor adding random spaces b
     * @param baseName the distribution name.
     * @param mean the mean for the distribution.
     * @param stddev the stddev for the distribution.
     * @param spacesDistr poisson distribution of random spaces to add in different parts of the descriptor.
     * @return never {@code null}.
     */
    private String composeDescriptionString(final String baseName, final double mean, final double stddev, PoissonDistribution spacesDistr) {
        return StringUtils.repeat(' ', spacesDistr.sample()) +
                baseName +
                StringUtils.repeat(' ', spacesDistr.sample()) +
                "(" +
                StringUtils.repeat(' ', spacesDistr.sample()) +
                mean +
                StringUtils.repeat(' ', spacesDistr.sample()) +
                ',' +
                StringUtils.repeat(' ', spacesDistr.sample()) +
                stddev +
                StringUtils.repeat(' ', spacesDistr.sample()) +
                ')' +
                StringUtils.repeat(' ', spacesDistr.sample());
    }

}
