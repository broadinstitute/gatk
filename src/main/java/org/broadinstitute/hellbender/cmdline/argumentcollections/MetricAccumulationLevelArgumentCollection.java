package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;

import java.io.Serializable;
import java.util.EnumSet;
import java.util.Set;

public class MetricAccumulationLevelArgumentCollection implements Serializable {

    private static final long serialVersionUID = 1L;

    @Argument(
            doc = "The level(s) at which to accumulate metrics. Possible values are {ALL_READS, SAMPLE, LIBRARY, READ GROUP}.",
            fullName = StandardArgumentDefinitions.METRIC_ACCUMULATION_LEVEL_LONG_NAME,
            shortName = StandardArgumentDefinitions.METRIC_ACCUMULATION_LEVEL_SHORT_NAME,
            optional = true
    )
    public Set<MetricAccumulationLevel> accumulationLevels = EnumSet.of(MetricAccumulationLevel.ALL_READS);
}
