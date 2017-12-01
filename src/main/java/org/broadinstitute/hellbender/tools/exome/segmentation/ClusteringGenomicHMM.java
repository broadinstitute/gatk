package org.broadinstitute.hellbender.tools.exome.segmentation;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.HMM;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Parent class for HMM-based segmentation of observed data type DATA and hidden data type HIDDEN.
 * For allele-fraction segmentation DATA is {@link AllelicCount}
 * and for copy-ratio segmentation DATA is Double (read counts).
 *
 * Hidden states are represented by their cluster indices, that is,
 * integers 0, 1, . . . K - 1 and store the values (of copy ratio or minor allele fraction) in a corresponding array.
 *
 * The model contains a memory length parameter representing the prior probability for the
 * state to be "forgotten" as a function of distance d between consecutive loci (targets or SNPs) --
 * the probability to remember a state is exp(-d/memoryLength).  If the state is forgotten, a new state
 * is chosen with probabilities given by an array of weights.
 *
 * Thus our transition probabilities are P(i -> j) = exp(-d/D) delta_{ij} + (1 - exp(-d/D)) weights[j]
 * where delta is the Kronecker delta and D is the memory length.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public abstract class ClusteringGenomicHMM<DATA, HIDDEN> implements HMM<DATA, SimpleInterval, Integer> {
    private final double memoryLength;
    private final List<HIDDEN> hiddenStateValues;
    private final int K;
    private final double logK;

    public ClusteringGenomicHMM(final List<HIDDEN> hiddenStateValues, final double memoryLength) {
        this.hiddenStateValues = new ArrayList<>(Utils.nonNull(hiddenStateValues));
        K = hiddenStateValues.size();
        logK = Math.log(K);
        Utils.nonEmpty(hiddenStateValues, "must have at least one hidden state");
        this.memoryLength = ParamUtils.isPositive(memoryLength, "CNV memory length must be positive");
    }

    // The following methods implement the HMM interface -------------------------------------------------
    @Override
    public List<Integer> hiddenStates() {
        return IntStream.range(0, hiddenStateValues.size()).boxed().collect(Collectors.toList());
    }

    @Override
    public double logPriorProbability(final Integer state, final SimpleInterval position) {
        return -logK;
    }

    // TODO: it's awkward that these are both required -- reason is that copy ratio emission is stored
    // TODO: by hidden state (Integer)
    // Child classes must specify their own emission likelihoods
    public abstract double logEmissionProbability(final DATA data, final Integer state, final SimpleInterval position);
    public abstract double logEmissionProbability(final DATA data, final HIDDEN hiddenStateValue);

    @Override
    public double logTransitionProbability(final Integer currentState, final SimpleInterval currentPosition,
                                           final Integer nextState, final SimpleInterval nextPosition) {
        return logTransitionProbability(currentState, nextState, calculateDistance(currentPosition, nextPosition));
    }

    @Override
    public void fillLogTransitionMatrix(final double[][] logTransitionMatrixBuffer, final SimpleInterval fromPosition, final SimpleInterval toPosition) {
        final double distance = calculateDistance(fromPosition, toPosition);
        final double decay = Math.exp(-distance / memoryLength);
        final double offDiagonalTerm = Math.log1p(-decay) - logK;
        final double diagonalTerm = Math.log(decay + (1 - decay)/K);
        for (int fromState = 0; fromState < K; fromState++) {
            for (int toState = 0; toState < K; toState++) {
                logTransitionMatrixBuffer[fromState][toState] = fromState == toState ? diagonalTerm : offDiagonalTerm;
            }
        }
    }
    // Done with implementation ----------------------------------------------------------------------------------------

    private double logTransitionProbability(final Integer currentState, final Integer nextState, final double distance) {
        final double pRemember = Math.exp(-distance / memoryLength);
        return Math.log((nextState.equals(currentState) ? pRemember : 0) + (1 - pRemember) / K);
    }

    public static double calculateDistance(final Locatable from, final Locatable to) {
        if (!from.getContig().equals(to.getContig())) {
            return Double.POSITIVE_INFINITY;
        } else {
            final double toMidpoint = (to.getStart() + to.getEnd())/2;
            final double fromMidpoint = (from.getStart() + from.getEnd())/2;
            return Math.abs(toMidpoint - fromMidpoint);
        }
    }

    public double getMemoryLength() { return memoryLength; }
    public HIDDEN getHiddenStateValue(final int state) { return hiddenStateValues.get(state); }
}
