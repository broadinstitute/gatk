package org.broadinstitute.hellbender.tools.exome.segmentation;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.HiddenMarkovModel;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Parent class for HMM-based segmentation of observed data type T.  For allele-fraction segmentation T is {@link AllelicCount}
 * and for copy-ratio segmentation T is Double (read counts).
 *
 * Hidden states are represented by their cluster indices, that is,
 * integers 0, 1, . . . K - 1 and store the values (of copy ratio or minor allele fraction) in a corresponding array.
 *
 * The model contains a memory length parameter representing the prior probability for the
 * state to be "forgotten" as a function of distance d between consecutive loci (targets or SNPs) --
 * the probability to remember a state is exp(-d/memoryLength).  If the state is forgotten, a new state
 * is chosen with probabilities given by an array of weights.
 *
 * Thus our transition probabilities are P(i -> j) = exp(-d/D) delta_{ij} + (1 - exp(-d/D) weights[j]
 * where delta is the Kronecker delta and D is the memory length.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public abstract class ClusteringGenomicHMM<T> implements HiddenMarkovModel<T, SimpleInterval, Integer> {
    private final double memoryLength;
    protected final double[] hiddenStateValues;
    protected final double[] weights;

    public ClusteringGenomicHMM(final double[] hiddenStateValues, final double[] weights, final double memoryLength) {
        Utils.nonNull(hiddenStateValues);
        Utils.nonNull(weights);
        Arrays.stream(weights).forEach(w -> ParamUtils.isPositiveOrZero(w, "weights may not be negative."));

        Utils.validateArg(hiddenStateValues.length == weights.length, "Must have one weight per minor allele fraction.");
        Utils.validateArg(hiddenStateValues.length > 0, "Must provide at least one minor allele fraction state.");

        this.hiddenStateValues = Arrays.copyOf(hiddenStateValues, hiddenStateValues.length);
        this.weights = MathUtils.normalizeFromRealSpace(weights);
        this.memoryLength = ParamUtils.isPositive(memoryLength, "CNV memory length must be positive");
    }

    // The following methods implement the HiddenMarkovModel interface -------------------------------------------------
    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, hiddenStateValues.length).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final SimpleInterval position) {
        return Math.log(weights[state]);
    }

    // TODO: it's awkward that these are both required -- reason is that copy ratio emission is stored
    // TODO: by hidden state (Integer)
    // Child classes must specify their own emission likelihoods
    @Override
    public abstract double logEmissionProbability(final T data, final Integer state, final SimpleInterval position);
    public abstract double logEmissionProbability(final T data, final double hiddenStateValue);

    @Override
    public double logTransitionProbability(final Integer currentState, final SimpleInterval currentPosition,
                                           final Integer nextState, final SimpleInterval nextPosition) {
        return logTransitionProbability(currentState, nextState, calculateDistance(currentPosition, nextPosition));
    }
    // Done with implementation ----------------------------------------------------------------------------------------

    private double logTransitionProbability(final Integer currentState, final Integer nextState, final double distance) {
        final double pRemember = Math.exp(-distance / memoryLength);
        return Math.log((nextState.equals(currentState) ? pRemember : 0) + (1 - pRemember)*weights[nextState]);
    }

    protected static double calculateDistance(final Locatable from, final Locatable to) {
        if (!from.getContig().equals(to.getContig())) {
            return Double.POSITIVE_INFINITY;
        } else {
            final double toMidpoint = (to.getStart() + to.getEnd())/2;
            final double fromMidpoint = (from.getStart() + from.getEnd())/2;
            return Math.abs(toMidpoint - fromMidpoint);
        }
    }

    public double getMemoryLength() { return memoryLength; }
    public double getWeight(final int k) { return weights[k]; }
    public double[] getWeights() { return Arrays.copyOf(weights, weights.length); }
    public double getHiddenStateValue(final int state) { return hiddenStateValues[state]; }
    public double[] getHiddenStateValues() { return Arrays.copyOf(hiddenStateValues, hiddenStateValues.length); }
}
