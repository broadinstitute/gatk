/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.downsampling;

import org.broadinstitute.gatk.utils.exceptions.UserException;

/**
 * Describes the method for downsampling reads at a given locus.
 */

public class DownsamplingMethod {
    /**
     * Type of downsampling to perform.
     */
    public final DownsampleType type;

    /**
     * Actual downsampling target is specified as an integer number of reads.
     */
    public final Integer toCoverage;

    /**
     * Actual downsampling target is specified as a fraction of total available reads.
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
     * Don't allow dcov values below this threshold for locus-based traversals (ie., Locus
     * and ActiveRegion walkers), as they can result in problematic downsampling artifacts
     */
    public static final int MINIMUM_SAFE_COVERAGE_TARGET_FOR_LOCUS_BASED_TRAVERSALS = 200;


    public DownsamplingMethod( DownsampleType type, Integer toCoverage, Double toFraction ) {
        this.type = type != null ? type : DEFAULT_DOWNSAMPLING_TYPE;

        if ( type == DownsampleType.NONE ) {
            this.toCoverage = null;
            this.toFraction = null;
        }
        else {
            this.toCoverage = toCoverage;
            this.toFraction = toFraction;
        }

        validate();
    }

    private void validate() {
        // Can't leave toFraction and toCoverage null unless type is NONE
        if ( type != DownsampleType.NONE && toFraction == null && toCoverage == null )
            throw new UserException("Must specify either toFraction or toCoverage when downsampling.");

        // Fraction and coverage cannot both be specified.
        if ( toFraction != null && toCoverage != null )
            throw new UserException("Downsampling coverage and fraction are both specified. Please choose only one.");

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
        StringBuilder builder = new StringBuilder("Downsampling Settings: ");

        if ( type == DownsampleType.NONE ) {
            builder.append("No downsampling");
        }
        else {
            builder.append(String.format("Method: %s, ", type));

            if ( toCoverage != null ) {
                builder.append(String.format("Target Coverage: %d", toCoverage));
            }
            else {
                builder.append(String.format("Target Fraction: %.2f", toFraction));
            }
        }

        return builder.toString();
    }
}
