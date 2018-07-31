package org.broadinstitute.hellbender.utils.read.markduplicates;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

/**
 * Abstract class that holds parameters and methods common to classes that optical duplicate detection.  We put them here so that
 * the explanation about how read names are parsed is in once place
 *
 * @author Tim Fennell
 */
public abstract class AbstractOpticalDuplicateFinderCommandLineProgram extends PicardCommandLineProgram {
    protected static Logger LOG = LogManager.getLogger(AbstractOpticalDuplicateFinderCommandLineProgram.class);

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    // The tool with which to find optical duplicates
    protected OpticalDuplicateFinder opticalDuplicateFinder = null;

    // Needed for testing
    public void setupOpticalDuplicateFinder() {
        this.opticalDuplicateFinder = new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX,
            opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE,null);//TODO logger firgure out, logger);
    }

    @Override
    protected String[] customCommandLineValidation() {
        setupOpticalDuplicateFinder();
        return super.customCommandLineValidation();
    }
}
