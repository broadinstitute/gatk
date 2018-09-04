package org.broadinstitute.hellbender.metrics;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.Serializable;

/**
 * MetricsArgumentCollection argument collection for QualityYield metrics. All members should be
 * instantiable as command line arguments.
 */
public class QualityYieldMetricsArgumentCollection extends MetricsArgumentCollection implements Serializable {

    public static final long serialVersionUID = 1L;

    @Argument(
            doc = "If available in the OQ tag, use the original quality scores " +
                    "as inputs instead of the quality scores in the QUAL field.",
            shortName = StandardArgumentDefinitions.USE_ORIGINAL_QUALITIES_SHORT_NAME,
            fullName = StandardArgumentDefinitions.USE_ORIGINAL_QUALITIES_LONG_NAME
    )
    public boolean useOriginalQualities = true;
}
