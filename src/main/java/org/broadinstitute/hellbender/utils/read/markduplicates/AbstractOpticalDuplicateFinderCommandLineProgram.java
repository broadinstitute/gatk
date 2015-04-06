package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;

/**
 * Abstract class that holds parameters and methods common to classes that optical duplicate detection.  We put them here so that
 * the explanation about how read names are parsed is in once place
 *
 * @author Tim Fennell
 */
public abstract class AbstractOpticalDuplicateFinderCommandLineProgram extends PicardCommandLineProgram {
    protected static Log LOG = Log.getInstance(AbstractOpticalDuplicateFinderCommandLineProgram.class);


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
            optional = true)
    public String READ_NAME_REGEX = OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX;

    @Argument(doc = "The maximum offset between two duplicte clusters in order to consider them optical duplicates. This " +
            "should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the " +
            "Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.")
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE;

    // The tool with which to find optical duplicates
    protected OpticalDuplicateFinder opticalDuplicateFinder = null;

    // Needed for testing
    public void setupOpticalDuplicateFinder() {
        this.opticalDuplicateFinder = new OpticalDuplicateFinder(READ_NAME_REGEX, OPTICAL_DUPLICATE_PIXEL_DISTANCE, LOG);
    }

    @Override
    protected String[] customCommandLineValidation() {
        setupOpticalDuplicateFinder();
        return super.customCommandLineValidation();
    }
}
