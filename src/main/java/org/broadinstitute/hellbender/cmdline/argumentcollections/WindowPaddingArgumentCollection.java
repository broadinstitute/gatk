package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;

/**
 * An abstract argument collection for window padding for slider walkers.
 *
 * Window padding is the length in each direction of the window to pad.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class WindowPaddingArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    /**
     * Get the window padding for making the windows
     */
    public abstract int getWindowPadding();
}
