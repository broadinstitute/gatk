package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.IntegerDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.LibraryStatistics;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Supported insert size distributions shapes.
 */
public enum InsertSizeDistributionShape {
    /**
     * The insert size distribution follow a normal distribution truncated at 0.
     */
    NORMAL("N", "Gauss", "Gaussian") {
        public AbstractIntegerDistribution fromMeanAndStdDeviation(final double mean, final double stddev) {
            final int seed = ((Double.hashCode(mean) * 31) + Double.hashCode(stddev) * 31 ) + this.name().hashCode() * 31;
            final RandomGenerator rdnGen = new JDKRandomGenerator();
            rdnGen.setSeed(seed);
            final RealDistribution normal = new NormalDistribution(mean, stddev);

            // We took the math for the truncted normal distribution from Wikipedia:
            // https://en.wikipedia.org/wiki/Truncated_normal_distribution
//            final NormalDistribution normal = new NormalDistribution(mean, stddev);
//            final double Z = 1.0 / (1.0 - normal.cumulativeProbability(-0.5));
//            expectedDensity = (x) -> Z * (normal.cumulativeProbability(x + 0.5) - normal.cumulativeProbability(x - 0.5));

            final double zeroCumulative = normal.cumulativeProbability(-0.5);
            final double zeroProb = normal.density(0);
            final double normalization = 1.0 / (1.0 - zeroCumulative); // = 1.0 / Z in Wikipedia
            final double mu = normal.getNumericalMean();
            final double sigmaSquare = normal.getNumericalVariance();
            final double sigma = Math.sqrt(sigmaSquare);
            final double newMean = mu + sigma * zeroProb * normalization;
            final double newVariance = sigmaSquare * (1 + normalization * (-mu / sigma) * zeroProb
                    - Math.pow(zeroProb * normalization, 2));
            // The actual distribution:
            return new AbstractIntegerDistribution(rdnGen) {

                private static final long serialVersionUID = -1L;

                @Override
                public double probability(int x) {
                    return  normalization * (normal.cumulativeProbability(x + 0.5) - normal.cumulativeProbability(x - 0.5));
                }

                @Override
                public double cumulativeProbability(int x) {
                    return (normal.cumulativeProbability(x + 0.5) - zeroCumulative) * normalization;
                }

                @Override
                public double getNumericalMean() {
                    return newMean;
                }

                @Override
                public double getNumericalVariance() {
                    return newVariance;
                }

                @Override
                public int getSupportLowerBound() {
                    return 0;
                }

                @Override
                public int getSupportUpperBound() {
                    return Integer.MAX_VALUE;
                }

                @Override
                public boolean isSupportConnected() {
                    return true;
                }
            };
        }

        @Override
        public AbstractIntegerDistribution fromSerializationFile(String whereFrom) {
            final ReadMetadata metaData = ReadMetadata.Serializer.readStandalone(whereFrom);
            final IntHistogram hist = new IntHistogram(2000);
            long modeCount = 0;
            for (final LibraryStatistics libStats : metaData.getAllLibraryStatistics().values()) {
                final IntHistogram.CDF cdf = libStats.getCDF();
                final long cdfTotalCount = cdf.getTotalObservations();
                final int size = cdf.size() - 1;
                hist.addObservations(0, Math.round(cdf.getFraction(0) * cdfTotalCount));
                for (int i = 1; i < size; i++) {
                    final double fraction = cdf.getFraction(i) - cdf.getFraction(i - 1);
                    final long count = Math.round(fraction * cdfTotalCount);
                    if (modeCount < count) {
                        modeCount = count;
                    }
                    hist.addObservations(i, count);
                }
            }
            if (hist.getTotalObservations() == 0) {
                throw new UserException.MalformedFile("Could not find any insert-sizes in " + whereFrom);
            }
            return hist.empiricalDistribution((int) Math.max(1, modeCount / 1_000_000));
        }


    }, /**
        * The insert size distribution follows a log-normal; i.e. the exp(isize) ~ Normal.
        */
      LOG_NORMAL("LogNormal", "LnN") {
        @Override
        public AbstractIntegerDistribution fromMeanAndStdDeviation(final double mean, final double stddev) {
            final double var = stddev * stddev;
            final double scale = 2 * Math.log(mean) - 0.5 * Math.log(var + mean * mean); // scale = mu in wikipedia article.
            final double shape = Math.sqrt(Math.log(1 + (var / (mean * mean)))); // shape = sigma in wikipedia article.
            final RealDistribution real =  new LogNormalDistribution(scale, shape);
            final int seed = (((Double.hashCode(mean) * 31) + Double.hashCode(stddev) * 31) + name().hashCode() * 31);
            final RandomGenerator rdnGen = new JDKRandomGenerator();
            rdnGen.setSeed(seed);
            return new AbstractIntegerDistribution(rdnGen) {

                private static final long serialVersionUID = -1L;

                @Override
                public double probability(int x) {
                    return real.cumulativeProbability(x + 0.5) - real.cumulativeProbability(x - 0.5);
                }

                @Override
                public double cumulativeProbability(int x) {
                    return real.cumulativeProbability(x + 0.5);
                }

                @Override
                public double getNumericalMean() {
                    return real.getNumericalMean();
                }

                @Override
                public double getNumericalVariance() {
                    return real.getNumericalVariance();
                }

                @Override
                public int getSupportLowerBound() {
                    return 0;
                }

                @Override
                public int getSupportUpperBound() {
                    return Integer.MAX_VALUE;
                }

                @Override
                public boolean isSupportConnected() {
                    return true;
                }
            };
        }
    }, /**
        * Arbitrary densities are set for each possible insert size.
        */
       EMPIRICAL("E", "Emp") {

        @Override
        public AbstractIntegerDistribution fromMeanAndStdDeviation(double mean, double stddev) {
            throw new UserException.BadInput("Empirical insert-size-distribution needs a meta-file");
        }

        @Override
        public AbstractIntegerDistribution fromTextFile(String whereFrom) {
            final IntHistogram hist = new IntHistogram(2000); // 2000 is the number of tracked values i.e. 0..2000
            long modeCount = 0;
            try (final BufferedReader reader = new BufferedReader(new InputStreamReader(BucketUtils.openFile(whereFrom)))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.startsWith(ReadMetadata.CDF_PREFIX)) {
                        final String[] cdf = line.substring(ReadMetadata.CDF_PREFIX.length() + 1).split("\t");
                        long leftCdf = 0;
                        for (int value = 0; value < cdf.length; value++) {
                            final long frequency = Long.parseLong(cdf[value]) - leftCdf;
                            hist.addObservations(value, frequency);
                            if (hist.getNObservations(value) > modeCount) {
                                modeCount = hist.getNObservations(value);
                            }
                        }
                    }
                }
                if (hist.getTotalObservations() == 0) {
                    throw new UserException.MalformedFile("Could not find any insert-sizes in " + whereFrom);
                }
                // We apply a smoothing that won't yield probabilities far below 10^-6 (Phred ~ 60)
                return hist.empiricalDistribution((int) Math.max(1, modeCount / 1_000_000));
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(whereFrom);
            }
        }
    };

    private final List<String> aliases;

    private static final Map<String, InsertSizeDistributionShape> byLowercaseName;

    static {
        final InsertSizeDistributionShape[] shapes = values();
        byLowercaseName = new HashMap<>(shapes.length * (1 + 5));
        for (final InsertSizeDistributionShape shape : shapes) {
            byLowercaseName.put(shape.name().toLowerCase(), shape);
            for (final String alias : shape.aliases) {
                byLowercaseName.put(alias.toLowerCase(), shape);
            }
        }
    }

    InsertSizeDistributionShape(final String ... names) {
        final List<String> nameList = Arrays.asList(names);
        aliases = Collections.unmodifiableList(nameList);
    }

    public List<String> aliases() {
        return aliases;
    }

    public static InsertSizeDistributionShape decode(final String name) {
        Utils.nonNull(name);
        return byLowercaseName.get(name.toLowerCase());
    }

    protected abstract AbstractIntegerDistribution fromMeanAndStdDeviation(final double mean, final double stddev);

    protected AbstractIntegerDistribution fromReadMetadataFile(final String whereFrom) {
        try {
            return fromSerializationFile(whereFrom);
        } catch (final RuntimeException ex) {
            return fromTextFile(whereFrom);
        }
    }

    protected AbstractIntegerDistribution fromTextFile(String whereFrom) {
        try (final BufferedReader reader = new BufferedReader(new InputStreamReader(BucketUtils.openFile(whereFrom)))) {
            String line;
            int value;
            double totalSum = 0;
            double totalSqSum = 0;
            long totalCount = 0;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(ReadMetadata.CDF_PREFIX)) {
                    final String[] cdf = line.substring(ReadMetadata.CDF_PREFIX.length() + 1).split("\t");
                    long leftCdf = 0;
                    for (value = 0; value < cdf.length; value++) {
                        final long frequency = Long.parseLong(cdf[value]) - leftCdf;
                        leftCdf += frequency;
                        totalSum += frequency * value;
                        totalSqSum += value * value * frequency;
                    }
                    totalCount += leftCdf;
                }
            }
            if (totalCount == 0) {
                throw new UserException.MalformedFile("Could not find any insert-sizes in " + whereFrom);
            }
            final double mean = totalSum / totalCount;
            final double stdDev = Math.sqrt(Math.abs(totalSqSum/totalCount - mean * mean));
            return fromMeanAndStdDeviation(mean, stdDev);
        } catch (final IOException ex2) {
            throw new UserException.CouldNotReadInputFile(whereFrom);
        } catch (final NumberFormatException ex2) {
            throw new UserException.MalformedFile("the CDF contains non-numbers in " + whereFrom);
        }
    }

    protected AbstractIntegerDistribution fromSerializationFile(String whereFrom) {
        final ReadMetadata metaData = ReadMetadata.Serializer.readStandalone(whereFrom);
        double totalSum = 0;
        double totalSqSum = 0;
        long totalCount = 0;
        for (final LibraryStatistics libStats : metaData.getAllLibraryStatistics().values()) {
            final IntHistogram.CDF cdf = libStats.getCDF();
            final long cdfTotalCount = cdf.getTotalObservations();
            final int size = cdf.size();
            for (int i = 1; i < size; i++) {
                final double fraction = cdf.getFraction(i) - cdf.getFraction(i - 1);
                final double count = fraction * cdfTotalCount;
                totalSum += count * i;
                totalSqSum += count * i * i;
            }
            totalCount += cdfTotalCount;
        }
        if (totalCount == 0) {
            throw new UserException.MalformedFile("Could not find any insert-sizes in " + whereFrom);
        }
        final double mean = totalSum / totalCount;
        final double variance = Math.abs(totalSqSum / totalCount - mean * mean);
        final double stdDev = Math.sqrt(variance);
        return fromMeanAndStdDeviation(mean, stdDev);
    }

}
