package org.broadinstitute.hellbender.cmdline.argumentcollections;

/**
 * Argument collection where the window padding is not exposed to the final user
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class DefaultWindowPaddingArgumentCollection extends WindowPaddingArgumentCollection {
    private static final long serialVersionUID = 1L;

    private final int windowPadding;

    public DefaultWindowPaddingArgumentCollection(int windowPadding) {
        if ( windowPadding < 0 ) {
            throw new IllegalArgumentException("Window padding must be >= 0");
        }
        this.windowPadding = windowPadding;
    }

    @Override
    public int getWindowPadding() {
        return windowPadding;
    }
}
