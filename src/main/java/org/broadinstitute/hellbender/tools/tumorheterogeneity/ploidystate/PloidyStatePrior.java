package org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PloidyStatePrior {
    private final int numPloidyStates;
    private final List<PloidyState> ploidyStates;
    private final Map<PloidyState, Double> logProbabilityMassFunctionMap;

    public PloidyStatePrior(final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap) {
        Utils.nonNull(unnormalizedLogProbabilityMassFunctionMap);
        Utils.validateArg(!unnormalizedLogProbabilityMassFunctionMap.isEmpty(), "Number of ploidy states must be positive.");
        numPloidyStates = unnormalizedLogProbabilityMassFunctionMap.size();
        ploidyStates = Collections.unmodifiableList(new ArrayList<>(unnormalizedLogProbabilityMassFunctionMap.keySet()));
        this.logProbabilityMassFunctionMap = normalize(new LinkedHashMap<>(unnormalizedLogProbabilityMassFunctionMap));
    }

    public int numPloidyStates() {
        return numPloidyStates;
    }

    public List<PloidyState> ploidyStates() {
        return ploidyStates;
    }

    public double logProbability(final PloidyState ploidyState) {
        Utils.nonNull(ploidyState);
        if (!logProbabilityMassFunctionMap.containsKey(ploidyState)) {
            throw new IllegalArgumentException("Ploidy-state prior not specified for given ploidy state.");
        }
        return logProbabilityMassFunctionMap.get(ploidyState);
    }

    //normalize log probabilities in mass function (which may be unnormalized)
    private Map<PloidyState, Double> normalize(final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap) {
        final List<PloidyState> states = new ArrayList<>(unnormalizedLogProbabilityMassFunctionMap.keySet());
        final double[] log10Probabilities = states.stream()
                .mapToDouble(s -> MathUtils.logToLog10(unnormalizedLogProbabilityMassFunctionMap.get(s))).toArray();
        final double[] probabilities = MathUtils.normalizeFromLog10ToLinearSpace(log10Probabilities);
        final LinkedHashMap<PloidyState, Double> logProbabilityMassFunctionMap = new LinkedHashMap<>();
        IntStream.range(0, unnormalizedLogProbabilityMassFunctionMap.size())
                .forEach(i -> logProbabilityMassFunctionMap.put(states.get(i), Math.log(probabilities[i])));
        return logProbabilityMassFunctionMap;
    }
}
