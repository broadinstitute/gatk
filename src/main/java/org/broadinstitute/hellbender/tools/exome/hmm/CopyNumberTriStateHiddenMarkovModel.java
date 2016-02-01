package org.broadinstitute.hellbender.tools.exome.hmm;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang.math.DoubleRange;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
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
 *         <dt>Mean event target distance</dt><dd>(D in the paper) the expectation of the distance (in base-pairs) between
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
     * Valid range for the event-rate parameter.
     */
    public static final DoubleRange VALID_EVENT_RATE_RANGE = new DoubleRange(0, .5);

    /**
     * Valid range for the mean number of targets per event parameter.
     */
    public static final DoubleRange VALID_MEAN_NUMBER_OF_TARGETS_PER_EVENT_RANGE = new DoubleRange(0, Double.MAX_VALUE);

    /**
     * Valid range for the deletion event distribution center parameter.
     */
    public static final DoubleRange VALID_DELETION_CENTER_RANGE = new DoubleRange(- Double.MAX_VALUE, 0);

    /**
     * Valid range for the duplication event distribution center parameter.
     */
    public static final DoubleRange VALID_DUPLICATION_CENTER_RANGE = new DoubleRange(0, Double.MAX_VALUE);

    /**
     * Valid range for the intra-event target distance.
     */
    public static final DoubleRange VALID_MEAN_INTRA_EVENT_TARGET_DISTANCE_RANGE
            = new DoubleRange(1, Double.MAX_VALUE);


    /**
     * The assumed stddev for coverage values in copy neutral segment.
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

    /**
     * Cached <code>log(p)</code> where <code>p</code> is the <i>event-rate</i>.
     */
    private final double logP;

    /**
     * Cached <code>log(1 - 2p)</code> where <code>p</code> is the <i>event-rate</i> .
     */
    private final double log1Minus2P;

    /**
     * Cached <code>log(q)</code> where <code>q</code> is the <i>inverse of the mean event size in targets</i>.
     */
    private final double logQ;

    /**
     * Cached <code>log(1 - q)</code> where <code>q</code> is the <i>inverse of the mean event size in targets</i>.
     */
    private final double log1MinusQ;

    /**
     * Cached <code>1/D</code> where <code>D</code> is the <i>average distance between targets in the same event</i>.
     */
    private final double invD;


    /**
     * Emission distribution by hidden state.
     */
    private final Map<CopyNumberTriState, NormalDistribution>
        emissionDistributionByState;

    /**
     * Holds the mean number of base-pairs between targets in the same event.
     */
    private final double D;

    /**
     * Holds the event-rate per target.
     */
    private final double P;

    /**
     * Mean number of targets per event.
     */
    private final double T;

    /**
     * Creates a new model instance.
     *
     * @param eventRate the event-rate for this instance.
     * @param meanNumberOfTargetsPerEvent the expectation for the number of targets in an event.
     * @param deletionMean the deletion depth of coverage negative shift.
     * @param duplicationMean the duplication depth of coverage positive shift.
     * @param meanEventTargetDistance the expectation of the distance between consecutive targets in an event.
     *
     * @throws IllegalArgumentException if any of the model parameters has an invalid value.
     */
    public CopyNumberTriStateHiddenMarkovModel(final double eventRate,
                                               final double meanNumberOfTargetsPerEvent,
                                               final double deletionMean,
                                               final double duplicationMean,
                                               final double meanEventTargetDistance) {
        ParamUtils.inRange(VALID_EVENT_RATE_RANGE, eventRate, "event rate");
        ParamUtils.inRange(VALID_MEAN_NUMBER_OF_TARGETS_PER_EVENT_RANGE, meanNumberOfTargetsPerEvent,
                "mean number of targets per event");
        ParamUtils.inRange(VALID_DELETION_CENTER_RANGE, deletionMean, "deletion mean depth coverage shift");
        ParamUtils.inRange(VALID_DUPLICATION_CENTER_RANGE, duplicationMean, "duplication mean depth coverage shift");
        ParamUtils.inRange(VALID_MEAN_INTRA_EVENT_TARGET_DISTANCE_RANGE, meanEventTargetDistance,
                "mean intra-event target distance");

        P = eventRate;
        logP = Math.log(eventRate);
        log1Minus2P = Math.log1p(-2 * eventRate);
        T = meanNumberOfTargetsPerEvent;
        logQ = -Math.log(meanNumberOfTargetsPerEvent);
        log1MinusQ = Math.log1p(-Math.exp(logQ));
        D = meanEventTargetDistance;
        invD = 1.0 / D;
        emissionDistributionByState = new EnumMap<>(CopyNumberTriState.class);
        emissionDistributionByState.put(CopyNumberTriState.NEUTRAL, new NormalDistribution(0, NEUTRAL_DISTRIBUTION_SD));
        emissionDistributionByState.put(CopyNumberTriState.DELETION, new NormalDistribution(deletionMean, DELETION_DISTRIBUTION_SD));
        emissionDistributionByState.put(CopyNumberTriState.DUPLICATION, new NormalDistribution(duplicationMean, DUPLICATION_DISTRIBUTION_SD));
    }

    @Override
    public List<CopyNumberTriState> hiddenStates() {
        return HIDDEN_STATES;
    }

    @Override
    public double logPriorProbability(final CopyNumberTriState state, final Target target) {
        Utils.nonNull(state);
        return state == CopyNumberTriState.NEUTRAL ? log1Minus2P : logP;
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
        if (currentState == CopyNumberTriState.NEUTRAL) {
            // inter-event or event start transition prob.:
            return nextState == CopyNumberTriState.NEUTRAL ? log1Minus2P : logP;
        } else {
            final double distance = calculateDistance(fromTarget, nextPosition);
            final double logF = -distance * invD;
            final double log1MinusF = Math.log1p(-Math.exp(logF));
            if (nextState == CopyNumberTriState.NEUTRAL) {
                // End of event transition prob.:
                return GATKProtectedMathUtils.naturalLogSumExp(logF + logQ, log1MinusF + log1Minus2P);
            } else if (currentState == nextState) {
                // Event extension transition prob.:
                return GATKProtectedMathUtils.naturalLogSumExp(logF + log1MinusQ, log1MinusF + logP);
            } else {
                // Event "sign switch" (del. <--> dup.) transition prob.:
                return log1MinusF + logP;
            }
        }
    }

    /**
     * Calculate the distance between two targets.
     * <p>
     * If any of the targets provided does not contain intervals, the distance is set to 0 in order to revert
     * to the non-distance dependent model described in Fromer et al. paper
     * <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3484655/table/tbl1/">Table 1</a>.
     * </p>
     * <p>
     * If both targets map to different chromosomes then we return {@link Double#POSITIVE_INFINITY}.
     * </p>
     * <p>
     * If both targets overlap, the distance returned is 0.0
     * </p>
     * <p>
     * If both targets are on the same contig and they don't overlap then the distance is the number of
     * bases between their two closest end points (either start or stop coordinates). This method
     * works regardless of the relative positioning of the from and to input targets.
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
            return 0.0;
        } else if (fromInterval.overlaps(toInterval)) {
            return 0.0;
        } if (!fromInterval.getContig().equals(toInterval.getContig())) {
            return Double.POSITIVE_INFINITY;
        } else if (fromInterval.getStart() < toInterval.getStart()) {
            return toInterval.getStart() - fromInterval.getEnd() - 1;
        } else {
            return fromInterval.getStart() - toInterval.getEnd() - 1;
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

    /**
     * Returns the event rate as the probability of starting an event
     * at any given target.
     *
     * @return a valid positive finite double value between 0 and .5
     */
    public double getEventRate() {
        return P;
    }

    /**
     * Average distance in bp between targets in the same event.
     * @return a valid finit positive value.
     */
    public double getMeanEventTargetDistance() {
        return D;
    }

    /**
     * Average number of targets per event.
     * @return a valid finite positive value.
     */
    public double getNumberOfTargetsInEvent() {
        return T;
    }
}
