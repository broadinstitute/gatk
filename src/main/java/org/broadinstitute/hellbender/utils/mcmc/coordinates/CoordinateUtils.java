package org.broadinstitute.hellbender.utils.mcmc.coordinates;

import org.apache.commons.math3.analysis.function.Logit;
import org.apache.commons.math3.analysis.function.Sigmoid;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Implements transformations between unbounded walker coordinates and bounded parameter values.
 * See Sec. 58.4 of the STAN manual at https://github.com/stan-dev/stan/releases/download/v2.12.0/stan-reference-2.12.0.pdf.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CoordinateUtils {
    private CoordinateUtils() {}

    public static double transformWalkerCoordinateToBoundedVariable(final double walkerCoordinate,
                                                                    final double min,
                                                                    final double max) {
        Utils.validateArg(min < max, "Minimum bound must be strictly less than maximum bound.");
        return min + (max - min) * new Sigmoid().value(walkerCoordinate);
    }

    public static double transformBoundedVariableToWalkerCoordinate(final double boundedVariable,
                                                                    final double min,
                                                                    final double max) {
        Utils.validateArg(min < max, "Minimum bound must be strictly less than maximum bound.");
        Utils.validateArg(min <= boundedVariable && boundedVariable <= max, "Variable value " + boundedVariable + " is not within variable bounds [" + min + ", " + max + "].");
        return new Logit().value((boundedVariable - min) / (max - min));
    }

    public static double calculateLogJacobianFactor(final double boundedVariable,
                                                    final double min,
                                                    final double max) {
        Utils.validateArg(min < max, "Minimum bound must be strictly less than maximum bound.");
        return boundedVariable < min || boundedVariable > max ?
                Double.NEGATIVE_INFINITY :
                FastMath.log(boundedVariable - min) + FastMath.log(max - boundedVariable) - FastMath.log(max - min);
    }
}
