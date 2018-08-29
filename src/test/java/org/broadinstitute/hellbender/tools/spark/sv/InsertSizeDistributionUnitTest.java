package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadataTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.IntToDoubleFunction;

/**
 * Unit tests for {@link InsertSizeDistribution}.
 */
public class InsertSizeDistributionUnitTest extends GATKBaseTest {

    private static final String READ_METADATA_FILE_NAME = "read-metadata-example.txt.gz";

    private static final File READ_METADATA_FILE = new File( publicTestDir + InsertSizeDistribution.class.getPackage().getName()
            .replace('.', File.separatorChar) + File.separatorChar + READ_METADATA_FILE_NAME);


    @Test
    public void testFromSerializedMetaData() throws IOException {
        final ReadMetadata readMetadata = ReadMetadataTest.composeTestReadMetadata();
        final File tempFile = File.createTempFile("test-rd", ".bin.gz");
        try {
            tempFile.delete();
            tempFile.deleteOnExit();
            ReadMetadata.Serializer.writeStandalone(readMetadata, tempFile.toString());
            final InsertSizeDistribution lnDist = new InsertSizeDistribution("LogNormal(" + tempFile.toString() + ")");
            Assert.assertEquals(lnDist.mean(), ReadMetadataTest.LIBRARY_STATISTICS_MEAN, ReadMetadataTest.LIBRARY_STATISTIC_MEAN_DIFF);
            Assert.assertEquals(lnDist.stddev(), ReadMetadataTest.LIBRARY_STATISTICS_SDEV, ReadMetadataTest.LIBRARY_STATISTICS_SDEV_DIFF);

        } finally {
            tempFile.delete();
        }
    }

    @Test
    public void testFromTextMetaData() throws IOException {
        final ReadMetadata readMetadata = ReadMetadataTest.composeTestReadMetadata();
        final File tempFile = File.createTempFile("test-rd", ".txt");
        try {
            tempFile.delete();
            tempFile.deleteOnExit();
            ReadMetadata.writeMetadata(readMetadata, tempFile.toString());
            final InsertSizeDistribution lnDist = new InsertSizeDistribution("LogNormal(" + tempFile.toString() + ")");
            Assert.assertEquals(lnDist.mean(), ReadMetadataTest.LIBRARY_STATISTICS_MEAN, ReadMetadataTest.LIBRARY_STATISTIC_MEAN_DIFF);
            Assert.assertEquals(lnDist.stddev(), ReadMetadataTest.LIBRARY_STATISTICS_SDEV, ReadMetadataTest.LIBRARY_STATISTICS_SDEV_DIFF);
        } finally {
            tempFile.delete();
        }
    }

    @Test
    public void testFromRealTextMetaData() {
        final InsertSizeDistribution lnDist  = new InsertSizeDistribution("LogNormal(" + READ_METADATA_FILE.toString() + ")");
        final InsertSizeDistribution nDist  = new InsertSizeDistribution("Normal(" + READ_METADATA_FILE.toString() + ")");
        for (final InsertSizeDistribution dist : Arrays.asList(lnDist, nDist)) {
            Assert.assertEquals(dist.mean(), 379.1432, 0.01); // calculated independently using R.
            Assert.assertEquals(dist.variance(), 18162., 2);
            Assert.assertEquals(dist.stddev(), 134.76, 0.02);
        }
    }

    @Test(dataProvider = "testData")
    public void testProbability(final String description, final int x, final double expected, final double logExpected) {
        final InsertSizeDistribution isd = new InsertSizeDistribution(description);
        Assert.assertEquals(isd.probability(x), expected, 0.00001);
        if (isd.probability(x) > 0 && expected > 0) {
            Assert.assertEquals(isd.logProbability(x), logExpected, 0.01);
        } else {
            //TODO currently we don't have the hability of producing a finite log-prob if the prob is == 0.
            //TODO due to limitations in apache common-math. So this else avoid to fail in these instances.
            Assert.assertTrue(Math.abs(isd.logProbability(x) - logExpected) < 0.01 || isd.logProbability(x) == Double.NEGATIVE_INFINITY || logExpected == Double.NEGATIVE_INFINITY);
        }
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
        for (final InsertSizeDistributionShape type : Arrays.asList(InsertSizeDistributionShape.NORMAL, InsertSizeDistributionShape.LOG_NORMAL)) {
            final List<String> distrNames = type.aliases();
            for (final double mean : means) {
                for (final double cv : cvs) {
                    final double stddev = mean * cv;
                    // We use alternative code to compose the densities using the formulas found in wikipedia articles.
                    // the actual implementation in main relies on apache common math, so is difficult to fall into the
                    // same error mode thus masking bugs.
                    final IntToDoubleFunction expectedDensity;
                    final IntToDoubleFunction expectedLogDensity;
                    if (type == InsertSizeDistributionShape.NORMAL) {
                        final NormalDistribution normal = new NormalDistribution(mean, stddev);
                        final double Z = 1.0 / (1.0 - normal.cumulativeProbability(-0.5));
                        expectedDensity = (x) -> Z * (normal.cumulativeProbability(x + 0.5) - normal.cumulativeProbability(x - 0.5));
                        expectedLogDensity = (x) -> Math.log(expectedDensity.applyAsDouble(x));
                    } else if (type == InsertSizeDistributionShape.LOG_NORMAL) {
                        final double var = stddev * stddev;
                        final double logMean = Math.log(mean) - Math.log(Math.sqrt(1 + (var / (mean * mean))));
                        final double logStddev = Math.sqrt(Math.log(1 + var / (mean * mean)));
                        final LogNormalDistribution logNormal = new LogNormalDistribution(logMean, logStddev);
                        expectedDensity = (x) -> logNormal.cumulativeProbability(x + 0.5) - logNormal.cumulativeProbability(x - 0.5);
                        expectedLogDensity = (x) -> Math.log(expectedDensity.applyAsDouble(x));
                    } else {
                        throw new IllegalStateException("test do not support one of the type supported by InsertSizeDistribution: " + type.aliases().get(0));
                    }
                    // We add fixed length cases
                    for (final int fixedSize : fixedSizes) {
                        final String distrName = distrNames.get(random.nextInt(distrNames.size()));
                        result.add(new Object[]{
                                composeDescriptionString(distrName, mean, stddev, spacesDistr),
                                fixedSize, expectedDensity.applyAsDouble(fixedSize),
                        expectedLogDensity.applyAsDouble(fixedSize)});
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
                        expectedLogDensity.applyAsDouble(x)});
                    }
                }
            }
        }
        // A couple of hard-wired cases to attest that the code above is generating genuine cases
        // rather than have the same error mode as the implementation in main.
        // I enclose commented out the code use to calculate these values in R:

        result.add(new Object[] { "N(300,150)", 231, 0.002447848, Math.log(0.002447848)});
        // > probs = pnorm(c(149.5, 151.5) , 300, 150)
        // > probs[2] - probs[1]
        // 0.002447848

        result.add(new Object[] { "lnN(100,10)", 103, 0.03655933, Math.log(0.03655933)});
        // > mu = 100; sigma = 10
        // > meanlog = log(mu/sqrt(1+ (sigma^2)/mu^2)) # eq. from wikipedia article.
        // > sdlog   = sqrt(exp(2*mu + sigma^2)*(exp(sigma^2)-1)) # eq. from wikipedia article.
        // > probs = plnorm(c(102.5, 103.5), meanlog, sdlog)
        // > probs[2] - probs[1]
        // 0.03655933

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
