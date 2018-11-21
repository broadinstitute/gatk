package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager;

import java.util.List;

/**
 * A basic interface for a class to be used with the StratificationManager system
 *
 * @author Mark DePristo
 * @since 3/28/12
 */
public interface Stratifier<T> {
    /**
     * @return a list of all objects states that may be provided by this States provider
     */
    public List<T> getAllStates();
}
