package org.broadinstitute.hellbender.tools.exome.hmm;

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

    public static final String EVENT_START_PROBABILITY_FULL_NAME = "eventStartProbability";
    public static final String EVENT_START_PROBABILITY_SHORT_NAME = "eventProb";
    public static final String MEAN_EVENT_SIZE_FULL_NAME = "meanEventSize";
    public static final String MEAN_EVENT_SIZE_SHORT_NAME = "eventSize";
    public static final String MEAN_DELETION_COVERAGE_SHIFT_FULL_NAME = "meanDeletionCoverageSigmaShift";
    public static final String MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME = "deletionShift";
    public static final String MEAN_DUPLICATION_COVERAGE_SHIFT_FULL_NAME = "meanDuplicationCoverageSigmaShift";
    public static final String MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME = "duplicationShift";

    @Argument(doc = "Probability that a base in a copy-neutral segment is followed by a base belonging to a CNV.",
              fullName = EVENT_START_PROBABILITY_SHORT_NAME,
              shortName = EVENT_START_PROBABILITY_FULL_NAME,
              optional = true)
    public double eventStartProbability = 1.0e-8;

    @Argument(doc = "Estimated mean size of non-copy-neutral segments, in base-pairs",
              fullName = MEAN_EVENT_SIZE_FULL_NAME,
              shortName = MEAN_EVENT_SIZE_SHORT_NAME,
              optional = true)
    public double meanEventSize = 70_000;

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
    public CopyNumberTriStateHiddenMarkovModelArgumentCollection() { }

    /**
     * Creates a new model instance using the current parameters.
     * @throws UserException.BadArgumentValue if any of the argument contain an invalid value.
     */
    public CopyNumberTriStateHiddenMarkovModel createModel() {
        return new CopyNumberTriStateHiddenMarkovModel(eventStartProbability, meanEventSize, meanDeletionCoverageShift,
                meanDuplicationCoverageShift);
    }
}
