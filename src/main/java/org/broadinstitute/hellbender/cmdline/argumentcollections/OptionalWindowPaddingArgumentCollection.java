package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Argument collection where the window padding is exposed to the final user as an optional argument
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class OptionalWindowPaddingArgumentCollection extends WindowPaddingArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName="windowPadding", shortName="windowPadding", doc = "", optional = true)
    protected int windowPadding;

    public OptionalWindowPaddingArgumentCollection(final int defaultWindowPadding) {
        if ( defaultWindowPadding < 0 ) {
            throw new IllegalArgumentException("Default window size must be >= 0");
        }
        windowPadding = defaultWindowPadding;
    }

    @Override
    public int getWindowPadding() {
        if ( windowPadding < 0 ) {
            throw new UserException.BadArgumentValue("--windowPadding", Integer.toString(windowPadding), "window padding must be >= 0");
        }
        return windowPadding;
    }
}
