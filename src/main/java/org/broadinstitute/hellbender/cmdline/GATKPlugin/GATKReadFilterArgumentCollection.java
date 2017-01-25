package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;

/**
 * An abstract ArgumentCollection for defining the set of read filter descriptor plugin arguments that are exposed to the user on the command line.
 *
 * Subclasses should provide {@link org.broadinstitute.barclay.argparser.Argument} annotations for the arguments that should be exposed.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @see org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor
 */
public abstract class GATKReadFilterArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Returns the enabled filter names. Order should be honored.
     *
     * Default implementation returns an empty list.
     */
    public Collection<String> getUserEnabledReadFilterNames() {
        return Collections.emptyList();
    }

    /**
     * Returns the disabled read filter names. Order should be honored.
     *
     * Default implementation returns an empty list.
     */
    public Collection<String> getUserDisabledReadFilterNames() {
        return Collections.emptyList();
    }

    /**
     * Returns {@code true} if all default filters are disabled; {@code false} otherwise.
     */
    public abstract boolean disableToolDefaultReadFilters();

}
