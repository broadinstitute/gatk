package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Argument collection where the window step is exposed to the final user as an optional argument
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class OptionalWindowStepArgumentCollection extends WindowStepArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName="windowStep", shortName="windowStep", doc = "", optional = true)
    protected int windowStep;

    /**
     * Initialize with a default value for the window size
     */
    public OptionalWindowStepArgumentCollection(final int defaultWindowStep) {
        if ( defaultWindowStep < 1 ) {
            throw new IllegalArgumentException("Default window step must be > 0");
        }
        windowStep = defaultWindowStep;
    }

    /**
     * @throws UserException.BadArgumentValue if the window size is smaller than 1
     */
    @Override
    public int getWindowStep() {
        if ( windowStep < 1 ) {
            throw new UserException.BadArgumentValue("--windowStep", Integer.toString(windowStep), "window size must be > 0");
        }
        return windowStep;
    }
}
