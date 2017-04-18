package org.broadinstitute.hellbender.utils.hmm;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A simple HMM model where the priors, transition and emission probabilities are all uninformative uniform
 * distributions.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class UninformativeTestHMM implements HMM<Integer, Integer, Integer> {

    private final int numStates;
    public UninformativeTestHMM(final int numStates) {
        ParamUtils.isPositive(numStates, "bad number of states");
        this.numStates = numStates;
    }

    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, numStates).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final Integer position) {
        Utils.validIndex(state, numStates);
        return - Math.log(numStates);
    }

    @Override
    public double logTransitionProbability(final Integer currentState, final Integer currentPosition, final Integer nextState,
                                           final Integer nextPosition) {
        Utils.validIndex(currentState, numStates);
        Utils.validIndex(nextState, numStates);
        return - Math.log(numStates);
    }

    @Override
    public double logEmissionProbability(final Integer data, final Integer state, final Integer position) {
        Utils.validIndex(state, numStates);
        return 0;
    }
}
