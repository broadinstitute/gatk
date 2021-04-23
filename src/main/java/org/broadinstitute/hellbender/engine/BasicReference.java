package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * A source of reference base calls.
 */
public interface BasicReference {
    /**
     * Get the reference base calls in the specified window.
     * The calls are encoded as ASCII characters, 1 per byte.
     * The are no guarantees about upper or lower case, or particular values that may or may not be present.
     */
    byte[] getBases( final SimpleInterval window );
}
