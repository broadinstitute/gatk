package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.tools.coveragemodel.XHMMTargetLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriStateTransitionProbabilityCache;
import org.broadinstitute.hellbender.utils.hmm.HiddenMarkovModel;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 *  This model has the following parameters:
 *  Event start probability: the probability of a CNV event beginning at any given base.
 *  Mean event size: the average size (in base-pairs) of an event.
 *  Mean coverage shift for deletions: the mean coverage in regions that have undergone a copy loss.
 *  Mean coverage shift for duplications:  the mean coverage in regions that have undergone a copy gain.
 *
 * <p>
 * This model assumes that target are provided in genomic coordinate order although it can handle either ascending
 * and descending orders.
 * </p>
 * <p>
 * When any of the targets provided does not have a defined interval, the model uses a default distance.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CopyNumberTriStateHiddenMarkovModel
        implements HiddenMarkovModel<Double, Target, CopyNumberTriState> {

    /**
     * Immutable list of hidden states.
     */
    private static final List<CopyNumberTriState> HIDDEN_STATES =
            Collections.unmodifiableList(Arrays.asList(CopyNumberTriState.values()));

    /**
    * The assumed stddev for likelihoods, regardless of copy number state.
    */
    private static final double EMISSION_SD = 1.0;

    private final XHMMTargetLikelihoodCalculator likelihoodsCalculator;

    private final CopyNumberTriStateTransitionProbabilityCache logTransitionProbabilityCache;

    //use this distance if target intervals are unspecified
    public static double DEFAULT_DISTANCE_BETWEEN_TARGETS = 10_000;

    //per-base probability of transition from neutral to CNV state
    private final double eventStartProbability;

    //average size in bases of CNVs
    private final double meanEventSize;

    private final static int RANDOM_SEED = 1767;

    /**
     * Creates a new model instance.
     *
     * @param eventStartProbability the probability per base pair of a transition from neutral to a CNV
     * @param meanEventSize the expectation of the distance between consecutive targets in an event
     * @param deletionMeanShift the deletion depth of coverage negative shift.
     * @param duplicationMeaShift the duplication depth of coverage positive shift.
     * @throws IllegalArgumentException if any of the model parameters has an invalid value.
     */
    public CopyNumberTriStateHiddenMarkovModel(final double eventStartProbability,
                                               final double meanEventSize, final double deletionMeanShift,
                                               final double duplicationMeaShift) {
        ParamUtils.inRange(eventStartProbability, 0, 1, "Event probability must be between 0 and 1.");
        ParamUtils.isNegativeOrZero(deletionMeanShift, "Deletion coverage shift must be negative.");
        ParamUtils.isPositiveOrZero(duplicationMeaShift, "Duplication coverage shift must be positive");
        ParamUtils.isPositive(meanEventSize, "Mean event size must be positive.");
        this.eventStartProbability = eventStartProbability;
        this.meanEventSize = meanEventSize;
        logTransitionProbabilityCache = new CopyNumberTriStateTransitionProbabilityCache(meanEventSize, eventStartProbability);
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
        likelihoodsCalculator = new XHMMTargetLikelihoodCalculator(deletionMeanShift, duplicationMeaShift, EMISSION_SD, rng);
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
     * Calculate the distance between two targets.
     * <p>
     * If any of the targets provided does not contain intervals, the distance is set to @{code DEFAULT_DISTANCE_BETWEEN_TARGETS}
     * in order to revert to a non-distance dependent model
     * </p>
     * <p>
     * If both targets map to different chromosomes then we return {@link Double#POSITIVE_INFINITY}.
     * </p>
     * <p>
     * Otherwise, the distance returned is the distance between their centers. This method
     * works regardless of the targets' relative positions.
     * </p>
     * @param fromTarget the previous target.
     * @param toTarget the next target.
     * @return any values between 0 and {@link Double#POSITIVE_INFINITY}.
     * @throws NullPointerException if any of the targets is {@code null}.
     */
    @VisibleForTesting
    double calculateDistance(final Target fromTarget, final Target toTarget) {
        final SimpleInterval fromInterval = fromTarget.getInterval();
        final SimpleInterval toInterval = toTarget.getInterval();
        if (fromInterval == null || toInterval == null) {
            return DEFAULT_DISTANCE_BETWEEN_TARGETS;
        } else if (!fromInterval.getContig().equals(toInterval.getContig())) {
            return Double.POSITIVE_INFINITY;
        } else {
            final double toMidpoint = (toInterval.getStart() + toInterval.getEnd())/2;
            final double fromMidpoint = (fromInterval.getStart() + fromInterval.getEnd())/2;
            return  Math.abs(toMidpoint - fromMidpoint);
        }
    }

    @Override
    public double logEmissionProbability(final Double data,
                                         final CopyNumberTriState state,
                                         final Target target) {
        return likelihoodsCalculator.logLikelihood(target, state.copyRatio, data);
    }

    @VisibleForTesting
    public Double randomDatum(final CopyNumberTriState state, final Random rdn) {
        Utils.nonNull(state);
        Utils.nonNull(rdn);
        return likelihoodsCalculator.generateRandomZScoreData(state.copyRatio);
    }

    /**
     * Returns the deletion mean shift.
     * @return a valid negative finite double value.
     */
    public double getDeletionMean() {
        return likelihoodsCalculator.deletionMean();
    }

    /**
     * Return the duplication mean shift.
     * @return a valid positive finite double value.
     */
    public double getDuplicationMean() {
        return likelihoodsCalculator.duplicationMean();
    }

    public double getEventStartProbability() { return eventStartProbability; }
    public double getMeanEventSize() { return meanEventSize; }
}
