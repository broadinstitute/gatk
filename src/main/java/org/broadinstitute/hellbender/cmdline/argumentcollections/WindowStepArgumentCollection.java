package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;

/**
 * An abstract argument collection for window step for slider walkers.
 *
 * Window step is the distance between the previous window and the next one.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class WindowStepArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    /**
     * Get the window step for making the windows
     */
    public abstract int getWindowStep();

}
