package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

/**
 * An interface for objects that know how to copy themselves deeply
 *
 * NOTE: Java already has a {@link Cloneable} interface, however, the great book of guidelines
 * for GATK4 developers says:
 *
 *      "Thou shalt not override clone() unless thou really know what thou is doing."
 *
 * The author is not sure about the latter.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public interface Duplicable {

    /**
     * Performs a deep copy
     *
     * @return a new instance of the object
     */
    Duplicable duplicate();

    /**
     * Whether or not the wrapped object is null
     *
     * @return boolean
     */
    boolean hasValue();

    /**
     * Returns the stored value
     */
    Object value();
}
