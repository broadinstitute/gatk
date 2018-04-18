package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import java.io.Serializable;
import java.util.List;

/**
 * An abstract ArgumentCollection for defining the set of annotation descriptor plugin arguments that are exposed to the user on the command line.
 *
 * Subclasses should provide {@link org.broadinstitute.barclay.argparser.Argument} annotations for the arguments that should be exposed.
 *
 * @see org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor
 */
public abstract class GATKAnnotationArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Returns the enabled annotation. Order should be honored.
     */
    public abstract List<String> getUserEnabledAnnotationNames();

    /**
     * Returns the enabled annotation. Order should be honored.
     */
    public abstract List<String> getUserEnabledAnnotationGroups();

    /**
     * Returns the disabled annotation names. Order should be honored.
     */
    public abstract List<String> getUserDisabledAnnotationNames();

    /**
     * Returns {@code true} if all default annotations are disabled; {@code false} otherwise.
     */
    public abstract boolean getDisableToolDefaultAnnotations();

    /**
     * Returns {@code true} if all annotations are enabled; {@code false} otherwise.
     */
    public abstract boolean getEnableAllAnnotations();
}