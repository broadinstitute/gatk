package org.broadinstitute.hellbender.utils.hmm;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Uniform emission model that takes the state priors and transitions as
 * {@link org.apache.commons.math3.linear.RealVector} and {@link RealMatrix}
 * instances respectively.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class FlatRealTestHMModel implements HiddenMarkovModel<Integer,Integer,Integer> {

    private final RealVector priors;
    private final RealMatrix transitions;
    private final int numStates;

    public FlatRealTestHMModel(final RealVector priors, final RealMatrix transitions) {
        this.priors = priors;
        this.transitions = transitions;
        Utils.validateArg(priors.getDimension() == this.transitions.getColumnDimension(), "bad dimensions");
        Utils.validateArg(priors.getDimension() == this.transitions.getRowDimension(), "bad dimensions");
        numStates = priors.getDimension();
    }

    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, numStates).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final Integer position) {
        return Math.log(priors.getEntry(state));
    }

    @Override
    public double logTransitionProbability(final Integer currentState,
                                           final Integer currentPosition,
                                           final Integer nextState,
                                           final Integer nextPosition) {
        return Math.log(transitions.getEntry(currentState, nextState));
    }

    @Override
    public double logEmissionProbability(final Integer data, final Integer state,
                                         final Integer position) {
        return 0;
    }
}
