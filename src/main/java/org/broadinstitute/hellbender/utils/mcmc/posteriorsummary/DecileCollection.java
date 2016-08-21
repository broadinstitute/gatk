package org.broadinstitute.hellbender.utils.mcmc.posteriorsummary;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Represents a set of deciles.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class DecileCollection implements Serializable {
    public enum ConstructionMode {
        SAMPLES, DECILES
    }

    private static final long serialVersionUID = 145L;

    public static final int NUM_DECILES = 11;
    private static final double EPSILON = 1E-10;

    private final Map<Decile, Double> deciles = new HashMap<>(NUM_DECILES);

    /**
     * Constructs a DecileCollection from a list of samples or deciles using Apache Commons {@link Percentile}.
     * @param values    list of samples or deciles (caution should be used if this contains NaN or infinite values)
     * @param mode      specifies whether values are samples or deciles
     */
    public DecileCollection(final List<Double> values, final ConstructionMode mode) {
        Utils.nonNull(values);
        Utils.validateArg(!values.isEmpty(), "Cannot construct deciles for empty list of samples.");

        final Decile[] decileKeys = Decile.values();

        if (mode == ConstructionMode.DECILES) {
            Utils.validateArg(values.size() == NUM_DECILES,
                    "DecileCollection construction mode was set to DECILES but an incorrect number of values was passed.");
            Utils.validateArg(Ordering.natural().isOrdered(values),
                    "DecileCollection construction mode was set to DECILES but an unsorted list of values was passed.");
            IntStream.range(0, NUM_DECILES).forEach(i -> deciles.put(decileKeys[i], values.get(i)));
        } else if (mode == ConstructionMode.SAMPLES) {
            final Percentile percentile = new Percentile();
            percentile.setData(Doubles.toArray(values));
            for (int i = 0; i < NUM_DECILES; i++) {
                //percentile.evaluate argument must be in (0, 100], so use EPSILON for 0th percentile
                final double decile = i == 0 ? percentile.evaluate(EPSILON) : percentile.evaluate(10 * i);
                deciles.put(decileKeys[i], decile);
            }
        }
    }

    /**
     * Gets the specified decile.
     */
    public double get(final Decile decile) {
        return deciles.get(decile);
    }

    /**
     * Gets a list of all deciles.
     */
    public List<Double> getAll() {
        return Collections.unmodifiableList(Stream.of(Decile.values()).map(deciles::get).collect(Collectors.toList()));
    }

    /**
     * Gets a list of deciles, excluding first and last.
     */
    public List<Double> getInner() {
        return Collections.unmodifiableList(Stream.of(Decile.values()).map(deciles::get).collect(Collectors.toList()).subList(1, NUM_DECILES - 1));
    }
}
