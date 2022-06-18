package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import java.io.Serializable;
import java.util.List;

/**
 * An abstract ArgumentCollection for defining the set of read transformer descriptor plugin arguments that are exposed to the user on the command line.
 *
 * Subclasses should provide {@link org.broadinstitute.barclay.argparser.Argument} annotations for the arguments that should be exposed.
 *
 * @author Dror Kessler (dror27)
 * @see GATKTransformerFilterPluginDescriptor
 */
public abstract class GATKReadTransformerArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Returns the enabled transformer names. Order should be honored.
     */
    public abstract List<String> getUserEnabledReadTransformerNames();

    /**
     * Returns the disabled read transformer names. Order should be honored.
     */
    public abstract List<String> getUserDisabledReadTransformerNames();

    /**
     * Returns {@code true} if all default transformers are disabled; {@code false} otherwise.
     */
    public abstract boolean getDisableToolDefaultReadTransformers();

}
