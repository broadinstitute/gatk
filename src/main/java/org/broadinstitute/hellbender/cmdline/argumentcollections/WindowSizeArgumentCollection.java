package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;

/**
 * An abstract argument collection for window size for slider walkers.
 *
 * Window size is the length of each window.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class WindowSizeArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    /**
     * Get the window size for making the windows
     */
    public abstract int getWindowSize();

}
