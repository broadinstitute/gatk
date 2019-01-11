package org.broadinstitute.hellbender.metrics;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.Serializable;

/**
 * Base class for defining a set of metrics collector arguments. The members should be
 * instantiable as command line arguments.
 */
public class MetricsArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "File to write the output to.")
    public String output;

}
