package org.broadinstitute.hellbender.utils.hmm;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A simple HMM model where the priors, transition and emission probabilities are all uninformative uniform
 * distributions.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class UninformativeTestHMModel implements HiddenMarkovModel<Integer, Integer, Integer> {

    private final int numStates;
    public UninformativeTestHMModel(final int numStates) {
        if (numStates <= 0) {
            throw new IllegalArgumentException("bad number of states");
        }
        this.numStates = numStates;
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
        return - Math.log(numStates);
    }

    @Override
    public double logTransitionProbability(final Integer currentState, final Integer currentPosition, final Integer nextState,
                                           final Integer nextPosition) {
        if (currentState < 0 || currentState >= numStates) {
            throw new IllegalArgumentException("bad state");
        } else if (nextState < 0 || nextState >= numStates) {
            throw new IllegalArgumentException("bad state");
        }
        return - Math.log(numStates);
    }

    @Override
    public double logEmissionProbability(final Integer data, final Integer state, final Integer position) {
        if (state < 0 || state >= numStates) {
            throw new IllegalArgumentException("bad state");
        }
        return 0;
    }
}
