package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.Serializable;

/**
 * Describes the method for downsampling reads at a given locus.
 */
public final class DownsamplingMethod implements Serializable {
    private static final long serialVersionUID = 1L;
    /**
     * Type of downsampling to perform.
     */
    public final DownsampleType type;

    /**
     * Actual downsampling target is specified as an integer number of reads.
     * This field is null if no downsampling is to be done,
     * or we are downsampling to a fraction.
     */
    public final Integer toCoverage;

    /**
     * Actual downsampling target is specified as a fraction of total available reads.
     * This field is null if no downsampling is to be done,
     * or we are downsampling to coverage.
     */
    public final Double toFraction;

    /**
     * Expresses no downsampling applied at all.
     */
    public static final DownsamplingMethod NONE = new DownsamplingMethod(DownsampleType.NONE, null, null);

    /**
     * Default type to use if no type is specified
     */
    public static final DownsampleType DEFAULT_DOWNSAMPLING_TYPE = DownsampleType.BY_SAMPLE;

    /**
     * If type is null, then DEFAULT_DOWNSAMPLING_TYPE will be used.
     * If type is not NONE, then either toFraction or toCoverage must be not null.
     * At most one of toFraction or toCoverage can be specified.
     * If specified, toCoverage must be > 0.
     * If specified, toFraction must be >= 0.0 and <= 1.0
     */
    public DownsamplingMethod(final DownsampleType type, final Integer toCoverage, final Double toFraction ) {
        this.type = type != null ? type : DEFAULT_DOWNSAMPLING_TYPE;

        if ( type == DownsampleType.NONE ) {
            this.toCoverage = null;
            this.toFraction = null;
        } else {
            this.toCoverage = toCoverage;
            this.toFraction = toFraction;
        }

        validate();
    }

    private void validate() {
        // Can't leave toFraction and toCoverage null unless type is NONE
        if ( type != DownsampleType.NONE && toFraction == null && toCoverage == null ) {
            throw new UserException("Must specify either toFraction or toCoverage when downsampling.");
        }

        // Fraction and coverage cannot both be specified.
        if ( toFraction != null && toCoverage != null ) {
            throw new UserException("Downsampling coverage and fraction are both specified. Please choose only one.");
        }

        // toCoverage must be > 0 when specified
        if ( toCoverage != null && toCoverage <= 0 ) {
            throw new UserException("toCoverage must be > 0 when downsampling to coverage");
        }

        // toFraction must be >= 0.0 and <= 1.0 when specified
        if ( toFraction != null && (toFraction < 0.0 || toFraction > 1.0) ) {
            throw new UserException("toFraction must be >= 0.0 and <= 1.0 when downsampling to a fraction of reads");
        }
    }

    public String toString() {
        final StringBuilder builder = new StringBuilder("Downsampling Settings: ");

        if ( type == DownsampleType.NONE ) {
            builder.append("No downsampling");
        } else {
            builder.append(String.format("Method: %s, ", type));

            if ( toCoverage != null ) {
                builder.append(String.format("Target Coverage: %d", toCoverage));
            } else {
                builder.append(String.format("Target Fraction: %.2f", toFraction));
            }
        }

        return builder.toString();
    }
}
