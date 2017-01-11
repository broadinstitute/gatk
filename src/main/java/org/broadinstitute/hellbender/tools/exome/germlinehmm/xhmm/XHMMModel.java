package org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.tools.coveragemodel.XHMMEmissionProbabilityCalculator;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriStateHiddenMarkovModel;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriStateTransitionProbabilityCache;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.HiddenMarkovModel;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 *  This XHMM model has the following parameters:
 *  <ul>
 *      <li>
 *          Event start probability: the probability of a CNV event beginning at any given base.
 *      </li>
 *      <li>
 *          Mean event size: the average size (in base-pairs) of an event.
 *      </li>
 *      <li>
 *          Mean coverage shift for deletions: the mean coverage in regions that have undergone a copy loss.
 *      </li>
 *      <li>
 *          Mean coverage shift for duplications:  the mean coverage in regions that have undergone a copy gain.
 *      </li>
 *  </ul>
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
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class XHMMModel extends CopyNumberTriStateHiddenMarkovModel<XHMMEmissionData> {

    private static final long serialVersionUID = 2702812961362183203L;

    /**
    * The assumed stddev for likelihoods, regardless of copy number state.
    */
    private static final double EMISSION_SD = 1.0;

    /**
     * A random seed for instantiating a random number generator for {@link XHMMEmissionProbabilityCalculator}
     */
    private static final int RANDOM_SEED = 1767;

    /**
     * Creates a new model instance.
     *
     * @param eventStartProbability the probability per base pair of a transition from neutral to a CNV
     * @param meanEventSize the expectation of the distance between consecutive targets in an event
     * @param deletionMeanShift the deletion depth of coverage negative shift.
     * @param duplicationMeaShift the duplication depth of coverage positive shift.
     * @throws IllegalArgumentException if any of the model parameters has an invalid value.
     */
    public XHMMModel(final double eventStartProbability,
                     final double meanEventSize,
                     final double deletionMeanShift,
                     final double duplicationMeaShift) {
        super(new XHMMEmissionProbabilityCalculator(
                deletionMeanShift,
                duplicationMeaShift,
                EMISSION_SD, RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED))),
                eventStartProbability,
                meanEventSize);
    }

    @VisibleForTesting
    public Double randomDatum(final CopyNumberTriState state, final Random rdn) {
        Utils.nonNull(state);
        Utils.nonNull(rdn);
        return ((XHMMEmissionProbabilityCalculator)targetLikelihoodCalculator).generateRandomZScoreData(state.copyRatio);
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

    /**
     * Returns the deletion mean shift.
     * @return a valid negative finite double value.
     */
    public double getDeletionMean() {
        return ((XHMMEmissionProbabilityCalculator) targetLikelihoodCalculator).deletionMean();
    }

    /**
     * Return the duplication mean shift.
     * @return a valid positive finite double value.
     */
    public double getDuplicationMean() {
        return ((XHMMEmissionProbabilityCalculator) targetLikelihoodCalculator).duplicationMean();
    }

    /**
     * Returns the event start probability.
     * @return a double between 0 and 1.
     */
    public double getEventStartProbability() { return eventStartProbability; }

    /**
     * Returns the mean size of events (in number of bases).
     * @return  a positive double.
     */
    public double getMeanEventSize() { return meanEventSize; }
}
