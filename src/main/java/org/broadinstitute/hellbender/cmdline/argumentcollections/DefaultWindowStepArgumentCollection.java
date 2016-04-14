package org.broadinstitute.hellbender.cmdline.argumentcollections;

/**
 * Argument collection where the window step is not exposed to the final user
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class DefaultWindowStepArgumentCollection extends WindowStepArgumentCollection {
    private static final long serialVersionUID = 1L;

    private final int windowStep;

    public DefaultWindowStepArgumentCollection(final int windowStep) {
        if ( windowStep < 1 ) {
            throw new IllegalArgumentException("Window step must be > 0");
        }
        this.windowStep = windowStep;
    }

    @Override
    public int getWindowStep() {
        return windowStep;
    }
}
