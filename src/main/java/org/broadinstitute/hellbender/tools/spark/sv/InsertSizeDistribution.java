package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.IntegerDistribution;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.Serializable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Holds the information characterizing and insert size distribution.
 */
public class InsertSizeDistribution implements Serializable {

    private static final long serialVersionUID = 1L;

    private static Pattern DESCRIPTION_PATTERN =
            Pattern.compile("^\\s*(?<name>[^\\s\\(\\)]+)\\s*\\((?<mean>[^,\\(\\)]+?)\\s*(?:,\\s*(?<stddev>[^,\\(\\)]+?)\\s*)?\\)\\s*");

    private final String description;

    private transient AbstractIntegerDistribution dist;

    private AbstractIntegerDistribution dist() {
        initializeDistribution();
        return dist;
    }

    public double mean() {
        return dist().getNumericalMean();
    }

    public double variance() { return dist().getNumericalVariance(); }

    public double stddev() { return Math.sqrt(variance()); }

    public InsertSizeDistribution(final String distrString) {
        this.description = distrString;
        initializeDistribution();
    }

    private void initializeDistribution() {
        if (!description.matches(DESCRIPTION_PATTERN.pattern())) {
            throw new UserException.BadInput("unsupported insert size distribution description format: " + description);
        }
        final Matcher matcher = DESCRIPTION_PATTERN.matcher(description);
        if (!matcher.find()) {
            throw new UserException.BadInput("the insert-size distribution spec is not up to standard: " + description);
        }
        final String nameString = matcher.group("name");
        final String meanString = matcher.group("mean");
        final String stddevString = matcher.group("stddev");
        final InsertSizeDistributionShape type = extractDistributionShape(nameString, description);
        if (stddevString != null) {
            final double mean = extractDoubleParameter("mean", description, meanString, 0, Double.MAX_VALUE);
            final double stddev = extractDoubleParameter("stddev", description, stddevString, 0, Double.MAX_VALUE);
            dist = type.fromMeanAndStdDeviation(mean, stddev);
        } else {
            dist = type.fromReadMetadataFile(meanString);
        }
    }

    private static InsertSizeDistributionShape extractDistributionShape(final String nameString, final String description) {
        final InsertSizeDistributionShape result = InsertSizeDistributionShape.decode(nameString);
        if (result == null) {
            throw new UserException.BadInput("unsupported insert size distribution name '" + nameString
                    + "' in description: " + description);
        }
        return result;
    }

    private static double extractDoubleParameter(final String name, final String description,
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
        if (!(obj instanceof InsertSizeDistribution)) {
            return false;
        } else {
            return ((InsertSizeDistribution)obj).dist().equals(dist());
        }
    }

    public int minimum() {
        return Math.max(0, dist().getSupportLowerBound());
    }

    public int maximum() {
        return Math.min(Integer.MAX_VALUE, dist().getSupportUpperBound());
    }

    public double probability(final int size) {
        return dist().probability(size);
    }

    public double logProbability(final int size) {
        return dist().logProbability(size);
    }
}
