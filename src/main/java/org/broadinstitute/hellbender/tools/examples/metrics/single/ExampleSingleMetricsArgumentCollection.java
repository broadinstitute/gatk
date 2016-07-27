package org.broadinstitute.hellbender.tools.examples.metrics.single;

import org.broadinstitute.hellbender.metrics.MetricsArgumentCollection;

import java.io.Serializable;

/**
 * Argument argument collection for Example single level metrics.
 */
public class ExampleSingleMetricsArgumentCollection extends MetricsArgumentCollection implements Serializable {

    public static final long serialVersionUID = 1L;

    // the example single level collector has no user-controlled arguments
}
