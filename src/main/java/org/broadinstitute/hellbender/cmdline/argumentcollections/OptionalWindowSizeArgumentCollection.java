package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Argument collection where the window size is exposed to the final user as an optional argument
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class OptionalWindowSizeArgumentCollection extends WindowSizeArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName="windowSize", shortName="windowSize", doc = "", optional = true)
    protected int windowSize;

    /**
     * Initialize with a default value for the window size
     */
    public OptionalWindowSizeArgumentCollection(final int defaultWindowSize) {
        if ( defaultWindowSize < 1 ) {
            throw new IllegalArgumentException("Default window size must be > 0");
        }
        windowSize = defaultWindowSize;
    }

    /**
     * @throws UserException.BadArgumentValue if the window size is smaller than 0
     */
    @Override
    public int getWindowSize() {
        if ( windowSize <= 0 ) {
            throw new UserException.BadArgumentValue("--windowSize", Integer.toString(windowSize), "window size must be > 0");
        }
        return windowSize;
    }
}
