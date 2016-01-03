package org.broadinstitute.hellbender.utils.hmm;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple test model where one state is "heavy" meaning that priors and transitions
 * to that state have large prob. that the rest of states (the light ones).
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HeavyStateTestHMModel implements HiddenMarkovModel<Integer, Integer, Integer> {

    private final int numStates;
    private final int heavyState;
    private final double heavyStateWeight;
    private final double lightStateProbs;

    public HeavyStateTestHMModel(final int numStates, final int heavyState, final double heavyStateWeight) {
        if (numStates <= 1) {
            throw new IllegalArgumentException("bad number of states");
        } else if (heavyState < 0 || heavyState >= numStates) {
            throw new IllegalArgumentException("bad heavy state");
        } else if (heavyStateWeight <= (1.0 / numStates) || heavyStateWeight >= 1.0) {
            throw new IllegalArgumentException("bad heavy state weight");
        }
        this.heavyState = heavyState;
        this.numStates = numStates;
        this.heavyStateWeight = heavyStateWeight;
        this.lightStateProbs = (1.0 - heavyStateWeight) / (double) (numStates - 1);
    }

    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, numStates).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final Integer position) {
        if (state < 0 || state >= numStates) {
            throw new IllegalArgumentException("bad state");
        }
        return heavyState == state ? Math.log(heavyStateWeight) : Math.log(lightStateProbs);
    }

    @Override
    public double logTransitionProbability(final Integer currentState, final Integer currentPosition, final Integer nextState,
                                           final Integer nextPosition) {
        if (currentState < 0 || currentState >= numStates) {
            throw new IllegalArgumentException("bad state");
        } else if (nextState < 0 || nextState >= numStates) {
            throw new IllegalArgumentException("bad state");
        }
        return nextState == heavyState ? Math.log(heavyStateWeight) : Math.log(lightStateProbs);
    }

    @Override
    public double logEmissionProbability(final Integer data, final Integer state, final Integer position) {
        if (state < 0 || state >= numStates) {
            throw new IllegalArgumentException("bad state");
        }
        return 0;
    }
}
