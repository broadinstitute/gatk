package org.broadinstitute.hellbender.tools.exome.hmm;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriStateTransitionProbabilityCache;
import org.broadinstitute.hellbender.utils.hmm.HiddenMarkovModel;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;

/**
 * Implements the Three state copy number distance dependent model described in
 * <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3484655">Fromer et al. Am J Hum Genet. 2012 Oct 5; 91(4):597-607</a>
 *
 * <p>
 *     This model has the following parameters:
 *     <dl>
 *         <dt>Event-Rate</dt><dd>(p in the paper) the exome-wide event rate; what is the probability of any given
 *         target to be the start target for an event (either a deletion or duplication).</dd>
 *         <dt>Mean number of targets per event</dt><dd>(T in the paper) the average number of targets per each event
 *         (either a deletion or duplication).</dd>
 *         <dt>Mean event target distance</dt><dd>(meanEventSize in the paper) the expectation of the distance (in base-pairs) between
 *         consecutive targets in an event.</dd>
 *         <dt>Mean coverage depth shift for deletions</dt><dd>(-M in the paper) what is the average negative shift
 *         in coverage in regions that have undergone a copy loss.</dd>
 *         <dt>Mean coverage depth shift for duplications</dt><dd>(+M in the paper) what is the average positive
 *         shift in coverage depth in regions that have undergone a copy gain</dd>
 *     </dl>
 *     <p>
 *      Notice that this implementation allows a different shift magnitude for deletions and duplications which is fix
 *      to the same value in the original paper.
 *     </p>
 * </p>
 *
 * <p>
 * This model assumes that target are provided in genomic coordinate order although it can handle either ascending
 * and descending orders.
 * </p>
 * <p>
 * When any of the targets provided does not have a defined interval, the model automatically switches to the
 * non-distance dependent simpler model described in the same publication.
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
    * The assumed stddev for coverage values in a copy neutral segment.
    */
    private static final double NEUTRAL_DISTRIBUTION_SD = 1.0;

    /**
     * The assumed stddev for coverage values in a deletion.
     */
    private static final double DELETION_DISTRIBUTION_SD = 1.0;

    /**
     * The assumed stddev for coverage values in a duplication.
     */
    private static final double DUPLICATION_DISTRIBUTION_SD = 1.0;

    private final Map<CopyNumberTriState, NormalDistribution> emissionDistributionByState;

    private final CopyNumberTriStateTransitionProbabilityCache logTransitionProbabilityCache;

    //use this distance if target intervals are unspecified
    public static double DEFAULT_DISTANCE_BETWEEN_TARGETS = 10_000;

    //per-base probability of transition from neutral to CNV state
    private final double eventStartProbability;

    //average size in bases of CNVs
    private final double meanEventSize;

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
        emissionDistributionByState = new EnumMap<>(CopyNumberTriState.class);
        emissionDistributionByState.put(CopyNumberTriState.NEUTRAL, new NormalDistribution(0, NEUTRAL_DISTRIBUTION_SD));
        emissionDistributionByState.put(CopyNumberTriState.DELETION, new NormalDistribution(deletionMeanShift, DELETION_DISTRIBUTION_SD));
        emissionDistributionByState.put(CopyNumberTriState.DUPLICATION, new NormalDistribution(duplicationMeaShift, DUPLICATION_DISTRIBUTION_SD));
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
        final NormalDistribution distribution = emissionDistributionByState.get(Utils.nonNull(state));
        return distribution.logDensity(data);
    }

    @VisibleForTesting
    public Double randomDatum(final CopyNumberTriState state, final Random rdn) {
        Utils.nonNull(state);
        Utils.nonNull(rdn);
        final NormalDistribution stateDistro = emissionDistributionByState.get(state);
        return rdn.nextGaussian() * stateDistro.getStandardDeviation() + stateDistro.getMean();
    }

    /**
     * Returns the deletion mean shift.
     * @return a valid negative finite double value.
     */
    public double getDeletionMean() {
        return emissionDistributionByState.get(CopyNumberTriState.DELETION).getMean();
    }

    /**
     * Return the duplication mean shift.
     * @return a valid positive finite double value.
     */
    public double getDuplicationMean() {
        return emissionDistributionByState.get(CopyNumberTriState.DUPLICATION).getMean();
    }

    public double getEventStartProbability() { return eventStartProbability; }
    public double getMeanEventSize() { return meanEventSize; }
}
