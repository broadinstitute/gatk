package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.EnumMap;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents a set of deciles.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class DecileCollection {
    private final EnumMap<Decile, Double> deciles = new EnumMap<>(Decile.class);

    /**
     * Constructs a DecileCollection from a list of samples using Apache Commons {@link Percentile}.
     * @param samples   list of samples (caution should be used if this contains NaN or infinite values)
     */
    public DecileCollection(final List<Double> samples) {
        Utils.nonNull(samples);
        Utils.validateArg(!samples.isEmpty(), "Cannot construct deciles for empty list of samples.");

        final Percentile percentile = new Percentile();
        percentile.setData(Doubles.toArray(samples));
        final Decile[] decileKeys = Decile.values();
        for (int i = 1; i < 10; i++) {
            final double decile = percentile.evaluate(10 * i);
            deciles.put(decileKeys[i - 1], decile);
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
}
