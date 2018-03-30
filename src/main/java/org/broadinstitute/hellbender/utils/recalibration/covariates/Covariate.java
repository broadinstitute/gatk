package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;

import java.io.Serializable;
import java.util.List;

/**
 * The Covariate interface. A Covariate is a feature used in the recalibration that can be picked out of the read.
 * In general most error checking and adjustments to the data are done before the call to the covariates getValue methods in order to speed up the code.
 * This unfortunately muddies the code, but most of these corrections can be done per read while the covariates get called per base, resulting in a big speed up.
 */
public interface Covariate extends Serializable {
    public static long serialVersionUID = 1L;

    /**
     * Initialize any member variables using the command-line arguments passed to the walker
     *
     * @param RAC the recalibration argument collection
     */
    public void initialize(final RecalibrationArgumentCollection RAC, final List<String> readGroups);

    /**
     * Calculates covariate values for all positions in the read.
     *
     * @param read   the read to calculate the covariates on.
     * @param header SAM header for the read
     * @param values the object to record the covariate values for every base in the read.
     * @param recordIndelValues indicates whether values of the covariate are to be recorded for indels
     */
    public void recordValues(final GATKRead read, final SAMFileHeader header, final PerReadCovariateMatrix values, final boolean recordIndelValues);

    /**
     * Converts the internal representation of the key to String format for file output.
     *
     * @param key the long representation of the key
     * @return a string representation of the key
     */
    public String formatKey(final int key);

    /**
     * Converts an Object key into a long key using only the lowest numberOfBits() bits
     *
     * Only necessary for on-the-fly recalibration when you have the object, but need to store it in memory in long format. For counting covariates
     * the getValues method already returns all values in long format.
     *
     * @param value the object corresponding to the covariate
     * @return a long representation of the object
     */
    public int keyFromValue(final Object value);

    /**
     * Returns the maximum value possible for any key representing this covariate
     *
     * @return the maximum value possible for any key representing this covariate
     */
    public int maximumKeyValue();

    /**
     * Returns the names of the covariate, which is the simple class name without the "Covariate" part;
     */
    default String parseNameForReport() {
        return getClass().getSimpleName().split("Covariate")[0];
    }

}

