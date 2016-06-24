package org.broadinstitute.hellbender.metrics;


import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;

/**
 * Base class for defining a set of metrics collector arguments. The members should be
 * instantiable as command line arguments.
 */
public class MetricsArgumentCollection {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to write the output to.")
    public File output;

}
