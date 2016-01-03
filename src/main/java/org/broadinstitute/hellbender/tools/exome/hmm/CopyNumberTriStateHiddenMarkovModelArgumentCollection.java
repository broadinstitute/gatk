package org.broadinstitute.hellbender.tools.exome.hmm;

import org.apache.commons.lang.math.DoubleRange;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * User arguments to compose a custom {@link CopyNumberTriStateHiddenMarkovModel}.
 *
 * <p>
 *     Default values are those that were show to work well with human data in paper
 *     <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3484655/">Fromer et al. 2012</a>.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CopyNumberTriStateHiddenMarkovModelArgumentCollection {

    public static final String EVENT_RATE_SHORT_NAME = "eventRate";
    public static final String EVENT_RATE_FULL_NAME = "eventSegmentOpenRate";
    public static final String MEAN_TARGETS_PER_SEGMENT_SHORT_NAME = "eventTargetCount";
    public static final String MEAN_TARGETS_PER_SEGMENT_FULL_NAME = "meanNumberOfTargetsPerSegment";
    public static final String MEAN_DISTANCE_BETWEEN_TARGETS_IN_SEGMENT_FULL_NAME = "meanDistanceBetweenSegmentTargets";
    public static final String MEAN_DISTANCE_BETWEEN_TARGETS_IN_SEGMENT_SHORT_NAME = "eventTargetDistance";
    public static final String MEAN_DELETION_COVERAGE_SHIFT_FULL_NAME = "meanDeletionCoverageSigmaShift";
    public static final String MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME = "deletionShift";
    public static final String MEAN_DUPLICATION_COVERAGE_SHIFT_FULL_NAME = "meanDuplicationCoverageSigmaShift";
    public static final String MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME = "duplicationShift";

    @Argument(doc = "Probability that a target is the first in a deletion or duplication segment",
              fullName = EVENT_RATE_FULL_NAME,
              shortName = EVENT_RATE_SHORT_NAME,
              optional = true)
    public double eventSegmentRate = 1.0e-8;

    @Argument(doc = "Estimated mean number of targets in non-copy-neutral segments",
              fullName = MEAN_TARGETS_PER_SEGMENT_FULL_NAME,
              shortName = MEAN_TARGETS_PER_SEGMENT_SHORT_NAME,
              optional = true)
    public double meanTargetsPerSegment = 6;

    @Argument(doc = "Estimated mean distance between contiguous targets in non-copy-neutral segments in base-pairs",
              fullName = MEAN_DISTANCE_BETWEEN_TARGETS_IN_SEGMENT_FULL_NAME,
              shortName = MEAN_DISTANCE_BETWEEN_TARGETS_IN_SEGMENT_SHORT_NAME,
              optional = true)
    public double meanDistanceBetweenTargetsInSegment = 70_000;

    @Argument(doc = "Estimated mean standardized coverage shift in deletion segments",
              fullName = MEAN_DELETION_COVERAGE_SHIFT_FULL_NAME,
              shortName = MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME,
              optional = true)
    public double meanDeletionCoverageShift = -3;

    @Argument(doc = "Estimated mean standardized coverage shift in duplication segments",
              fullName = MEAN_DUPLICATION_COVERAGE_SHIFT_FULL_NAME,
              shortName = MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME,
              optional = true)
    public double meanDuplicationCoverageShift = 3;

    /**
     * Creates a new model argument collection taking on the default values suggested in
     *  <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3484655/">Fromer et al. 2012</a>.
     */
    public CopyNumberTriStateHiddenMarkovModelArgumentCollection() {

    }

    public void validate() {
        checkParameterRange(CopyNumberTriStateHiddenMarkovModel.VALID_EVENT_RATE_RANGE, eventSegmentRate, EVENT_RATE_FULL_NAME);
        checkParameterRange(CopyNumberTriStateHiddenMarkovModel.VALID_DELETION_CENTER_RANGE, meanDeletionCoverageShift, MEAN_DELETION_COVERAGE_SHIFT_FULL_NAME);
        checkParameterRange(CopyNumberTriStateHiddenMarkovModel.VALID_DUPLICATION_CENTER_RANGE, meanDuplicationCoverageShift, MEAN_DUPLICATION_COVERAGE_SHIFT_FULL_NAME);
        checkParameterRange(CopyNumberTriStateHiddenMarkovModel.VALID_MEAN_INTRA_EVENT_TARGET_DISTANCE_RANGE, meanDistanceBetweenTargetsInSegment, MEAN_DISTANCE_BETWEEN_TARGETS_IN_SEGMENT_FULL_NAME);
        checkParameterRange(CopyNumberTriStateHiddenMarkovModel.VALID_MEAN_NUMBER_OF_TARGETS_PER_EVENT_RANGE, meanTargetsPerSegment, MEAN_TARGETS_PER_SEGMENT_FULL_NAME);
    }

    private void checkParameterRange(final DoubleRange validRange, final double value, final String argumentName) {
        if (!validRange.containsDouble(value)) {
            throw new UserException.BadArgumentValue(argumentName, String.valueOf(value),
                    String.format("must be in the range [%g, %g]",
                    validRange.getMinimumDouble(), validRange.getMaximumDouble()));
        }
    }

    /**
     * Creates a new model instance using the current parameters.
     * @throws UserException.BadArgumentValue if any of the argument contain an invalid value.
     */
    public CopyNumberTriStateHiddenMarkovModel createModel() {
        validate();
        return new CopyNumberTriStateHiddenMarkovModel(eventSegmentRate, meanTargetsPerSegment, meanDeletionCoverageShift,
                meanDuplicationCoverageShift, meanDistanceBetweenTargetsInSegment);
    }
}
