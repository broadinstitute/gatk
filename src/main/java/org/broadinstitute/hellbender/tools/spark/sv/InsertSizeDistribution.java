package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Holds the information characterizing and insert size distribution.
 */
public class InsertSizeDistribution implements Serializable {

    private static final long serialVersionUID = -1L;

    @VisibleForTesting
    static final Type[] SUPPORTED_TYPES = { new NormalType(), new LogNormalType() };

    public interface Type {
        List<String> getNames();

        AbstractRealDistribution fromMeanAndStdDeviation(final double mean, final double stddev);
    }

    public static class NormalType implements Type {

        @Override
        public List<String> getNames() {
            return Collections.unmodifiableList(Arrays.asList("Normal", "N", "Norm", "Gauss", "Gaussian"));
        }

        @Override
        public AbstractRealDistribution fromMeanAndStdDeviation(final double mean, final double stddev) {
            return new NormalDistribution(mean, stddev);
        }
    }

    public static class LogNormalType implements Type {

        @Override
        public List<String> getNames() {
            return Collections.unmodifiableList(Arrays.asList("logN", "lnN", "logNorm", "lnNorm", "logNormal", "lnNormal"));
        }

        @Override
        public AbstractRealDistribution fromMeanAndStdDeviation(final double mean, final double stddev) {
            final double var = stddev * stddev;
            final double scale = 2 * Math.log(mean) - 0.5 * Math.log(var + mean * mean); // scale = mu in wikipedia article.
            final double shape = Math.sqrt(Math.log(1 + (var / (mean * mean)))); // shape = sigma in wikipedia article.
            return new LogNormalDistribution(scale, shape);
        }
    }

    private static Pattern DESCRIPTION_PATTERN =
            Pattern.compile("^\\s*(?<name>[^\\s\\(\\)]+)\\s*\\((?<mean>[^,\\(\\)]+?)\\s*,\\s*(?<stddev>[^,\\(\\)]+?)\\s*\\)\\s*");

    private final String description;

    private transient AbstractRealDistribution dist;

    private AbstractRealDistribution dist() {
        initializeDistribution();
        return dist;
    }

    public double average() {
        return dist().getNumericalMean();
    }

    public InsertSizeDistribution(final String distrString) {
        this.description = distrString;
        initializeDistribution();
    }

    private void initializeDistribution() {
        if (!description.matches(DESCRIPTION_PATTERN.pattern())) {
            throw new UserException.BadInput("unsupported insert size distribution description format: " + description);
        }
        final Matcher matcher = DESCRIPTION_PATTERN.matcher(description);
        matcher.find();
        final String nameString = matcher.group("name");
        final String meanString = matcher.group("mean");
        final String stddevString = matcher.group("stddev");
        final Type type = extractDistributionType(nameString, description);
        final double mean = extractDoubleParameter("mean", description, meanString, 0, Double.MAX_VALUE);
        final double stddev = extractDoubleParameter("stddev", description, stddevString, 0, Double.MAX_VALUE);
        dist = type.fromMeanAndStdDeviation(mean, stddev);
    }

    private static final Type extractDistributionType(final String nameString, final String description) {
        for (final Type candidate : SUPPORTED_TYPES) {
            if (candidate.getNames().stream().anyMatch(name -> name.toLowerCase().equals(nameString.trim().toLowerCase()))) {
                return candidate;
            }
        }
        throw new UserException.BadInput("unsupported insert size distribution name '" + nameString
                + "' in description: " + description);
    }

    private static final double extractDoubleParameter(final String name, final String description,
                                                       final String valueString, final double min,
                                                       final double max) {
        final double value;
        try {
            value = Double.parseDouble(valueString);
        } catch (final NumberFormatException ex) {
            throw new UserException.BadInput("bad " + name + " string '" + valueString + "' in insert size distribution description: "  + description);
        }
        if (value < min || Double.isInfinite(value) || Double.isNaN(value) || value > max) {
            throw new UserException.BadInput("bad " + name + " value '" + value + "' in insert size distribution description: " + description);
        }
        return value;
    }

    @Override
    public String toString() {
        return description;
    }

    @Override
    public int hashCode() {
        return dist().hashCode();
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == null || !(obj instanceof InsertSizeDistribution)) {
            return false;
        } else {
            return ((InsertSizeDistribution)obj).dist().equals(dist());
        }
    }

    public int minimum() {
        return (int) Math.max(0, dist().getSupportLowerBound());
    }

    public int maximum() {
        return (int) Math.min(Integer.MAX_VALUE, dist().getSupportUpperBound());
    }

    public double density(final int size) {
        return dist().density(size);
    }

    public double logDensity(final int size) {
        return dist().logDensity(size);
    }
}
