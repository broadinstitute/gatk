package org.broadinstitute.hellbender.utils.mcmc.coordinates;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Represents a walker position in unbounded N-dimensional space.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class WalkerPosition extends ArrayList<Double> {
    private static final long serialVersionUID = 79452354L;
    private final int numDimensions;

    public WalkerPosition(final List<Double> walkerPosition) {
        super(Collections.unmodifiableList(Utils.nonNull(walkerPosition)));
        Utils.validateArg(!walkerPosition.isEmpty(), "Dimension of walker space must be greater than zero.");
        numDimensions = walkerPosition.size();
    }

    public int numDimensions() {
        return numDimensions;
    }
}
