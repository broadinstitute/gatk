package org.broadinstitute.hellbender.tools.exome.eval;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Collection of user arguments to tune evaluation filters.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class EvaluationFiltersArgumentCollection implements Cloneable {

    public static final String MINIMUM_TRUTH_SEGMENT_LENGTH_SHORT_NAME = "minTruthLen";
    public static final String MINIMUM_TRUTH_SEGMENT_LENGTH_FULL_NAME = "minimumTruthSegmentLength";
    public static final String MINIMUM_CALLED_SEGMENT_LENGTH_SHORT_NAME = "minCallLen";
    public static final String MINIMUM_CALLED_SEGMENT_LENGTH_FULL_NAME = "minimumCalledSegmentLength";
    public static final String MINIMUM_TRUTH_SEGMENT_QUALITY_SHORT_NAME = "minTruthQual";
    public static final String MINIMUM_TRUTH_SEGMENT_QUALITY_FULL_NAME = "minimumTruthSegmentQuality";
    public static final String MINIMUM_CALLED_SEGMENT_QUALITY_SHORT_NAME = "minCallQual";
    public static final String MINIMUM_CALLED_SEGMENT_QUALITY_FULL_NAME = "minimumCalledSegmentQuality";
    public static final String APPLY_MULTI_ALLELIC_TRUTH_FILTER_SHORT_NAME = "applyMATFilter";
    public static final String APPLY_MULTI_ALLELIC_TRUTH_FILTER_FULL_NAME = "applyMultiAllelicTruthFilter";
    public static final String APPLY_MULTI_ALLELIC_CALLED_FILTER_SHORT_NAME = "applyMACFilter";
    public static final String APPLY_MULTI_ALLELIC_CALLED_FILTER_FULL_NAME = "applyMultiAllelicCalledFilter";
    public static final String MAXIMUM_TRUTH_EVENT_FREQUENCY_SHORT_NAME = "maxTruthFreq";
    public static final String MAXIMUM_TRUTH_EVENT_FREQUENCY_FULL_NAME = "maximumTruthEventFrequency";
    public static final String MAXIMUM_CALLED_EVENT_FREQUENCY_SHORT_NAME = "maxCallFreq";
    public static final String MAXIMUM_CALLED_EVENT_FREQUENCY_FULL_NAME = "maximumCalledEventFrequency";

    public static final int DEFAULT_MINIMUM_TRUTH_SEGMENT_LENGTH = 1;
    public static final int DEFAULT_MINIMUM_CALLED_SEGMENT_LENGTH = 1;
    public static final double DEFAULT_MINIMUM_TRUTH_SEGMENT_QUALITY = 0.0;
    public static final double DEFAULT_MINIMUM_CALLED_SEGMENT_QUALITY = 0.0;
    public static final double DEFAULT_MAXIMUM_TRUTH_EVENT_FREQUENCY = 1.0;
    public static final double DEFAULT_MAXIMUM_CALLED_EVENT_FREQUENCY = 1.0;

    public EvaluationFiltersArgumentCollection() {
    }

    @Argument(
            doc = "Minimum Truth segment length to consider the segment callable.",
            shortName = MINIMUM_TRUTH_SEGMENT_LENGTH_SHORT_NAME,
            fullName = MINIMUM_TRUTH_SEGMENT_LENGTH_FULL_NAME,
            optional = true
    )
    public int minimumTruthSegmentLength = DEFAULT_MINIMUM_TRUTH_SEGMENT_LENGTH;

    @Argument(
            doc = "Minimum called segment length to consider the segment trust-worthy.",
            shortName = MINIMUM_CALLED_SEGMENT_LENGTH_SHORT_NAME,
            fullName = MINIMUM_CALLED_SEGMENT_LENGTH_FULL_NAME,
            optional = true
    )
    public int minimumCalledSegmentLength = DEFAULT_MINIMUM_CALLED_SEGMENT_LENGTH;

    @Argument(
            doc = "Minimum Truth segment quality.",
            shortName = MINIMUM_TRUTH_SEGMENT_QUALITY_SHORT_NAME,
            fullName = MINIMUM_TRUTH_SEGMENT_QUALITY_FULL_NAME,
            optional = true
    )
    public double minimumTruthSegmentQuality = DEFAULT_MINIMUM_TRUTH_SEGMENT_QUALITY;

    @Argument(
            doc = "Minimum called segment quality.",
            shortName = MINIMUM_CALLED_SEGMENT_QUALITY_SHORT_NAME,
            fullName = MINIMUM_CALLED_SEGMENT_QUALITY_FULL_NAME,
            optional = true
    )
    public double minimumCalledSegmentQuality = DEFAULT_MINIMUM_CALLED_SEGMENT_QUALITY;

    @Argument(
            doc = "Apply the multi-allelic truth segment filter",
            shortName = APPLY_MULTI_ALLELIC_TRUTH_FILTER_SHORT_NAME,
            fullName = APPLY_MULTI_ALLELIC_TRUTH_FILTER_FULL_NAME,
            optional = true
    )
    public boolean applyMultiAllelicTruthFilter = false;

    @Argument(
            doc = "Apply the multi-allelic called segment filter",
            shortName = APPLY_MULTI_ALLELIC_CALLED_FILTER_SHORT_NAME,
            fullName = APPLY_MULTI_ALLELIC_CALLED_FILTER_FULL_NAME,
            optional = true
    )
    public boolean applyMultiAllelicCalledFilter = false;

    @Argument(
            doc = "Maximum frequency for truth events at a position",
            shortName = MAXIMUM_TRUTH_EVENT_FREQUENCY_SHORT_NAME,
            fullName = MAXIMUM_TRUTH_EVENT_FREQUENCY_FULL_NAME,
            optional = true
    )
    public double maximumTruthEventFrequency = DEFAULT_MAXIMUM_TRUTH_EVENT_FREQUENCY;

    @Argument(
            doc = "Maximum frequency for truth events at a position",
            shortName = MAXIMUM_CALLED_EVENT_FREQUENCY_SHORT_NAME,
            fullName = MAXIMUM_CALLED_EVENT_FREQUENCY_FULL_NAME,
            optional = true
    )
    public double maximumCalledEventFrequency = DEFAULT_MAXIMUM_CALLED_EVENT_FREQUENCY;

    @Override
    public EvaluationFiltersArgumentCollection clone() {
        try {
            return (EvaluationFiltersArgumentCollection) super.clone();
        } catch (final CloneNotSupportedException e) {
            throw new GATKException.ShouldNeverReachHereException("this class is cloneable", e);
        }
    }
}
