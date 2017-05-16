package org.broadinstitute.hellbender.utils.hmm;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple test model where one state is "heavy" meaning that priors and transitions
 * to that state have large prob. that the rest of states (the light ones).
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HeavyStateTestHMM implements HMM<Integer, Integer, Integer> {

    private final int numStates;
    private final int heavyState;
    private final double heavyStateWeight;
    private final double lightStateProbs;

    public HeavyStateTestHMM(final int numStates, final int heavyState, final double heavyStateWeight) {
        Utils.validateArg(numStates > 1, "bad number of states");
        Utils.validIndex(heavyState, numStates);
        ParamUtils.inRange(heavyStateWeight, (1.0 / numStates), 1.0, "bad heavy state weight");
        this.heavyState = heavyState;
        this.numStates = numStates;
        this.heavyStateWeight = heavyStateWeight;
        lightStateProbs = (1.0 - heavyStateWeight) / (double) (numStates - 1);
    }

    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, numStates).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final Integer position) {
        Utils.validIndex(state, numStates);
        return heavyState == state ? Math.log(heavyStateWeight) : Math.log(lightStateProbs);
    }

    @Override
    public double logTransitionProbability(final Integer currentState, final Integer currentPosition, final Integer nextState,
                                           final Integer nextPosition) {
        Utils.validIndex(currentState, numStates);
        Utils.validIndex(nextState, numStates);
        return nextState == heavyState ? Math.log(heavyStateWeight) : Math.log(lightStateProbs);
    }

    @Override
    public double logEmissionProbability(final Integer data, final Integer state, final Integer position) {
        Utils.validIndex(state, numStates);
        return 0;
    }
}
