package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.Serializable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by valentin on 5/18/17.
 */
public class InsertSizeDistribution implements Serializable {

    private static final long serialVersionUID = 1L;

    public static final String MEAN_DISTR_PARAM_NAME = "mean";
    public static final String SD_DISTR_PARAM_NAME = "sd";

    private static Pattern NORMAL_DISTR_PATTERN = Pattern.compile(
            String.format("^N(?:orm(?:al)?)?\\(\\s*(?<%s>\\S+)\\s*,\\s*(?<%s>\\S+)\\s*\\)",MEAN_DISTR_PARAM_NAME, SD_DISTR_PARAM_NAME));

    private final String description;

    private transient RealDistribution dist;

    private RealDistribution dist() {
        initializeDistribution();
        return dist;
    }

    public InsertSizeDistribution(final String distrString) {
        this.description = distrString;
        initializeDistribution();
    }

    private void initializeDistribution() {
        if (dist == null) {
            if (description.matches(NORMAL_DISTR_PATTERN.pattern())) {
                dist = parseNormalDistribution(description);
            } else {
                throw new IllegalArgumentException("unsupported insert size distribution description: " + description);
            }
        }
    }

    private static RealDistribution parseNormalDistribution(final String distrString) {
        final Matcher matcher = NORMAL_DISTR_PATTERN.matcher(distrString);
        if (!matcher.find()) {
            throw new IllegalArgumentException("bad description: " + distrString);
        }
        final double mean = Double.parseDouble(matcher.group(MEAN_DISTR_PARAM_NAME));
        final double sd = Double.parseDouble(matcher.group(SD_DISTR_PARAM_NAME));
        if (!Double.isFinite(mean)) {
            throw new UserException.BadInput(String.format("bad insert distribution mean value, must be a finite double but you provided: %d",  mean));
        } else if (!Double.isFinite(sd) || sd <= 0) {
            throw new UserException.BadInput(String.format("bad insert distribution std. dev. value must be a strictly positive finite but you provided: %d", sd));
        } else {
            return new NormalDistribution(mean, sd);
        }
    }

    public int minimum() {
        return (int) Math.max(0, dist().getSupportLowerBound());
    }

    public int maximum() {
        return (int) Math.min(Integer.MAX_VALUE, dist().getSupportUpperBound());
    }

    public double probability(final int size) {
        return dist().density(size);
    }

    public double logProbability(final int size) {
        return Math.log(probability(size));
    }
}
