package org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate;

import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents the prior on a set of {@link PloidyState} objects.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PloidyStatePrior {
    private final List<PloidyState> ploidyStates;
    private final Map<PloidyState, Double> ploidyStateToLogProbabilityMap;
    private final Map<Integer, Double> copyNumberToLogProbabilityMap;
    private final int maxCopyNumber;

    public PloidyStatePrior(final Map<PloidyState, Double> ploidyStateToUnnormalizedLogProbabilityMap) {
        Utils.nonNull(ploidyStateToUnnormalizedLogProbabilityMap);
        Utils.validateArg(!ploidyStateToUnnormalizedLogProbabilityMap.isEmpty(), "Number of ploidy states must be positive.");
        ploidyStates = Collections.unmodifiableList(new ArrayList<>(ploidyStateToUnnormalizedLogProbabilityMap.keySet()));
        ploidyStateToLogProbabilityMap = normalize(new LinkedHashMap<>(ploidyStateToUnnormalizedLogProbabilityMap));
        copyNumberToLogProbabilityMap = calculateCopyNumberToLogProbabilityMap(ploidyStateToLogProbabilityMap);
        maxCopyNumber = ploidyStates.stream().mapToInt(PloidyState::total).max().getAsInt();
    }

    public List<PloidyState> ploidyStates() {
        return ploidyStates;
    }

    public int maxCopyNumber() {
        return maxCopyNumber;
    }

    public double logProbability(final PloidyState ploidyState) {
        Utils.nonNull(ploidyState);
        if (!ploidyStateToLogProbabilityMap.containsKey(ploidyState)) {
            throw new IllegalArgumentException("Prior not specified for given ploidy state.");
        }
        return ploidyStateToLogProbabilityMap.get(ploidyState);
    }

    public double logProbability(final int copyNumber) {
        if (!copyNumberToLogProbabilityMap.containsKey(copyNumber)) {
            throw new IllegalArgumentException("Prior not specified for given copy number.");
        }
        return copyNumberToLogProbabilityMap.get(copyNumber);
    }

    //marginalize over ploidy states with the same total copy number
    private static Map<Integer, Double> calculateCopyNumberToLogProbabilityMap(final Map<PloidyState, Double> ploidyStateToLogProbabilityMap) {
        final Map<Integer, List<Double>> copyNumberToLogProbabilitiesMap = new HashMap<>();
        for (final PloidyState ploidyState : ploidyStateToLogProbabilityMap.keySet()) {
            final int copyNumber = ploidyState.total();
            if (!copyNumberToLogProbabilitiesMap.keySet().contains(copyNumber)) {
                copyNumberToLogProbabilitiesMap.put(copyNumber, new ArrayList<>());
            }
            final double logProbability = ploidyStateToLogProbabilityMap.get(ploidyState);
            copyNumberToLogProbabilitiesMap.get(copyNumber).add(logProbability);
        }
        return normalize(copyNumberToLogProbabilitiesMap.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, e -> GATKProtectedMathUtils.logSumExp(e.getValue()))));
    }

    //normalize log probabilities in map
    private static <T> Map<T, Double> normalize(final Map<T, Double> stateToUnnormalizedLogProbabilityMap) {
        final List<T> states = new ArrayList<>(stateToUnnormalizedLogProbabilityMap.keySet());
        final double[] log10Probabilities = states.stream()
                .mapToDouble(s -> MathUtils.logToLog10(stateToUnnormalizedLogProbabilityMap.get(s))).toArray();
        final double[] probabilities = MathUtils.normalizeFromLog10ToLinearSpace(log10Probabilities);
        return IntStream.range(0, stateToUnnormalizedLogProbabilityMap.size()).boxed()
                .collect(Collectors.toMap(states::get, i -> Math.log(probabilities[i])));
    }
}
