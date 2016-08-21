package org.broadinstitute.hellbender.utils.mcmc.coordinates;

import org.apache.commons.math3.analysis.function.Logit;
import org.apache.commons.math3.analysis.function.Sigmoid;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a position constrained to lie on the unit simplex and implements transformations to and from
 * unbounded walker coordinates.
 * See Sec. 58.6 of the STAN manual at https://github.com/stan-dev/stan/releases/download/v2.12.0/stan-reference-2.12.0.pdf.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class SimplexPosition extends ArrayList<Double> {
    private static final long serialVersionUID = 79454L;
    private static final double NORMALIZATION_EPSILON = 1E-6;
    private static final double EPSILON = 1E-10;
    private final int numDimensions;

    public SimplexPosition(final List<Double> simplexPosition) {
        super(Collections.unmodifiableList(new ArrayList<>(Utils.nonNull(simplexPosition))));
        Utils.validateArg(simplexPosition.size() > 1, "Number of dimensions must be strictly greater than 1.");
        final double normalization = simplexPosition.stream().mapToDouble(Double::doubleValue).sum();
        Utils.validateArg(Math.abs(1. - normalization) <= NORMALIZATION_EPSILON,
                "Position does not lie on the simplex within a tolerance of NORMALIZATION_EPSILON.");
        numDimensions = simplexPosition.size();
    }
    
    public int numDimensions() {
        return numDimensions;
    }

    public static double calculateLogJacobianFactor(final SimplexPosition simplexPosition) {
        final List<Double> breakProportions = calculateBreakProportionsFromSimplexPosition(simplexPosition);
        return IntStream.range(0, simplexPosition.size() - 1).boxed()
                .mapToDouble(i -> FastMath.log(Math.max(EPSILON, simplexPosition.get(i))) + FastMath.log(Math.max(EPSILON, 1. - breakProportions.get(i)))).sum();
    }

    public static WalkerPosition calculateWalkerPositionFromSimplexPosition(final SimplexPosition simplexPosition) {
        final List<Double> breakProportions = calculateBreakProportionsFromSimplexPosition(simplexPosition);
        return new WalkerPosition(calculateWalkerPositionFromBreakProportions(breakProportions));
    }

    public static SimplexPosition calculateSimplexPositionFromWalkerPosition(final WalkerPosition walkerPosition) {
        final List<Double> breakProportions = calculateBreakProportionsFromWalkerPosition(walkerPosition);
        final List<Double> simplexPosition = calculateSimplexPositionFromBreakProportions(breakProportions);
        return new SimplexPosition(simplexPosition);
    }

    private static List<Double> calculateSimplexPositionFromBreakProportions(final List<Double> breakProportions) {
        final int numDimensions = breakProportions.size() + 1;
        final List<Double> simplexPosition = new ArrayList<>();
        double cumulativeSum = 0.;
        for (int dimensionIndex = 0; dimensionIndex < numDimensions - 1; dimensionIndex++) {
            final double fraction = (1. - cumulativeSum) * breakProportions.get(dimensionIndex);
            simplexPosition.add(fraction);
            cumulativeSum += fraction;
        }
        simplexPosition.add(1. - cumulativeSum);
        return new SimplexPosition(simplexPosition);
    }

    private static List<Double> calculateBreakProportionsFromSimplexPosition(final SimplexPosition simplexPosition) {
        final int numDimensions = simplexPosition.numDimensions();
        final List<Double> breakProportions = new ArrayList<>();
        double cumulativeSum = 0.;
        for (int dimensionIndex = 0; dimensionIndex < numDimensions - 1; dimensionIndex++) {
            final double breakProportion = simplexPosition.get(dimensionIndex) / (1. - cumulativeSum);
            breakProportions.add(breakProportion);
            cumulativeSum += simplexPosition.get(dimensionIndex);
        }
        return breakProportions;
    }

    private static List<Double> calculateBreakProportionsFromWalkerPosition(final WalkerPosition walkerPosition) {
        final int numDimensions = walkerPosition.size() + 1;
        return IntStream.range(0, numDimensions - 1).boxed()
                .map(i -> new Sigmoid().value(walkerPosition.get(i) + FastMath.log(1. / (numDimensions - (i + 1)))))
                .collect(Collectors.toList());
    }

    private static List<Double> calculateWalkerPositionFromBreakProportions(final List<Double> breakProportions) {
        final int numDimensions = breakProportions.size() + 1;
        return IntStream.range(0, numDimensions - 1).boxed()
                .map(i -> new Logit().value(breakProportions.get(i)) - FastMath.log(1. / (numDimensions - (i + 1))))
                .map(x -> Math.max(-Double.MAX_VALUE, Math.min(Double.MAX_VALUE, x)))
                .collect(Collectors.toList());
    }
}
