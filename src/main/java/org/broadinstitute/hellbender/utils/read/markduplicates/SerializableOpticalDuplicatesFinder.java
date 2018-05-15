package org.broadinstitute.hellbender.utils.read.markduplicates;

import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import java.io.Serializable;

/**
 * An intermediate class flagged as being serializable so that the Picard OpticalDuplicateFinder can be serialized despite
 * it not being marked as such. Since there is no internal state beyond initialization settings in the picard OpticalDuplicates
 * finder this should not cause serialization problems.
 */
public class SerializableOpticalDuplicatesFinder extends OpticalDuplicateFinder implements Serializable{
    private static final long serialVersionUID = 1L;

    public SerializableOpticalDuplicatesFinder(String read_name_regex, int optical_duplicate_pixel_distance) {
        super(read_name_regex, optical_duplicate_pixel_distance, null);
    }

    public SerializableOpticalDuplicatesFinder() {
        super();
    }
}
