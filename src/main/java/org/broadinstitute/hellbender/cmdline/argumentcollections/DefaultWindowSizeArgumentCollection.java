package org.broadinstitute.hellbender.cmdline.argumentcollections;

/**
 * Argument collection where the window size is not exposed to the final user
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class DefaultWindowSizeArgumentCollection extends WindowSizeArgumentCollection {
    private static final long serialVersionUID = 1L;

    private final int windowSize;

    public DefaultWindowSizeArgumentCollection(final int windowSize) {
        if ( windowSize < 1 ) {
            throw new IllegalArgumentException("Window size must be > 0");
        }
        this.windowSize = windowSize;
    }

    @Override
    public int getWindowSize() {
        return windowSize;
    }
}
