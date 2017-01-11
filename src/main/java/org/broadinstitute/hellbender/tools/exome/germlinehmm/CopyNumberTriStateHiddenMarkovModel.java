package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.tools.coveragemodel.TargetLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.coveragemodel.XHMMEmissionProbabilityCalculator;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.HiddenMarkovModel;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 *  A generic copy number tri-state HMM has the following parameters:
 *  <ul>
 *      <li>
 *          Emission probability calculator: an implementation of {@link TargetLikelihoodCalculator}
 *      </li>
 *      <li>
 *          Event start probability: the probability of a CNV event beginning at any given base.
 *      </li>
 *      <li>
 *          Mean event size: the average size (in base-pairs) of an event.
 *      </li>
 *  </ul>
 *
 * <p>
 * This model assumes that target are provided in genomic coordinate order although it can handle either ascending
 * and descending orders.
 * </p>
 * <p>
 * When any of the targets provided does not have a defined interval, the model uses a default distance
 * {@link #DEFAULT_DISTANCE_BETWEEN_TARGETS}
 * </p>
 *
 * @param <D> is the emission data type
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class CopyNumberTriStateHiddenMarkovModel<D>
        implements HiddenMarkovModel<D, Target, CopyNumberTriState>, Serializable {

    private static final long serialVersionUID = -3367028114020534823L;

    /**
     * Immutable list of hidden states.
     */
    protected static final List<CopyNumberTriState> HIDDEN_STATES =
            Collections.unmodifiableList(Arrays.asList(CopyNumberTriState.values()));

    protected final TargetLikelihoodCalculator<D> targetLikelihoodCalculator;

    protected final CopyNumberTriStateTransitionProbabilityCache logTransitionProbabilityCache;

    /* use this distance if target intervals are unspecified */
    public static double DEFAULT_DISTANCE_BETWEEN_TARGETS = 10_000;

    /* per-base probability of transition from neutral to CNV state */
    protected final double eventStartProbability;

    /* average size in bases of CNVs */
    protected final double meanEventSize;

    /**
     * Creates a new model instance.
     *
     * @param emissionProbabilityCalculator emission probability calculator for each target
     * @param eventStartProbability the probability per base pair of a transition from neutral to a CNV
     * @param meanEventSize the expectation of the distance between consecutive targets in an event
     * @throws IllegalArgumentException if any of the model parameters has an invalid value.
     */
    public CopyNumberTriStateHiddenMarkovModel(@Nonnull final TargetLikelihoodCalculator<D> emissionProbabilityCalculator,
                                               final double eventStartProbability,
                                               final double meanEventSize) {
        ParamUtils.inRange(eventStartProbability, 0, 1, "Event probability must be between 0 and 1.");
        ParamUtils.isPositive(meanEventSize, "Mean event size must be positive.");
        this.eventStartProbability = eventStartProbability;
        this.meanEventSize = meanEventSize;
        logTransitionProbabilityCache = new CopyNumberTriStateTransitionProbabilityCache(meanEventSize, eventStartProbability);
        this.targetLikelihoodCalculator = Utils.nonNull(emissionProbabilityCalculator,
                "The emission probability calculator must be non-null");
    }

    @Override
    public List<CopyNumberTriState> hiddenStates() {
        return HIDDEN_STATES;
    }

    @Override
    public double logPriorProbability(final CopyNumberTriState state, final Target target) {
        Utils.nonNull(state);
        return logTransitionProbabilityCache.logProbability(Integer.MAX_VALUE, state, CopyNumberTriState.NEUTRAL);
    }

    /**
     * {@inheritDoc}
     *
     * @param currentState the destination state.
     * @param fromTarget the target before the transition.
     * @param nextState the source state.
     * @param nextPosition the target after the transition.
     * @return {@inheritDoc}.
     * @throws IllegalArgumentException if any of the inputs (target or state) is {@code null}.
     */
    @Override
    public double logTransitionProbability(final CopyNumberTriState currentState, final Target fromTarget,
                                           final CopyNumberTriState nextState, final Target nextPosition) {
        Utils.nonNull(currentState);
        Utils.nonNull(nextState);
        Utils.nonNull(fromTarget);
        Utils.nonNull(nextPosition);
        final double distance = calculateDistance(fromTarget, nextPosition);
        return logTransitionProbabilityCache.logProbability((int) distance, nextState, currentState);
    }

    /**
     * See {@link Target#calculateDistance(Target, Target, double)}
     *
     * @param fromTarget first target
     * @param toTarget second target
     * @return distance
     */
    public static double calculateDistance(final Target fromTarget, final Target toTarget) {
        return Target.calculateDistance(fromTarget, toTarget, DEFAULT_DISTANCE_BETWEEN_TARGETS);
    }

    @Override
    public double logEmissionProbability(final D emissionData,
                                         final CopyNumberTriState state,
                                         final Target target) {
        return targetLikelihoodCalculator.logLikelihood(emissionData, state.copyRatio, target);
    }

    public double getEventStartProbability() { return eventStartProbability; }

    public double getMeanEventSize() { return meanEventSize; }
}
