package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.io.Serializable;


/**
 * An argument collection for use with tools that mark optical
 * duplicates.
 */
public final class OpticalDuplicatesArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String READ_NAME_REGEX_LONG_NAME = "read-name-regex";
    public static final String OPTICAL_DUPLICATE_PIXEL_DISTANCE_LONG_NAME = "optical-duplicate-pixel-distance";

    @Argument(doc = "Regular expression that can be used to parse read names in the incoming SAM file. Read names are " +
             "parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used " +
             "to estimate the rate of optical duplication in order to give a more accurate estimated library size. " +
             "Set this option to null to disable optical duplicate detection. " +
             "The regular expression should contain three capture groups for the three variables, in order. " +
             "It must match the entire read name. " +
             "Note that if the default regex is specified, a regex match is not actually done, but instead the read name " +
             " is split on colon character. " +
             "For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. " +
             "For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values.",
             fullName = READ_NAME_REGEX_LONG_NAME,
             optional = true)
    public String READ_NAME_REGEX = OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX;

    @Argument(doc = "The maximum offset between two duplicate clusters in order to consider them optical duplicates. This " +
             "should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the " +
             "Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.",
              fullName = OPTICAL_DUPLICATE_PIXEL_DISTANCE_LONG_NAME,
              optional = true)
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE;
}
