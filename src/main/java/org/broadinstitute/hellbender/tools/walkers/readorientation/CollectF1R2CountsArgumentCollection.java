package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.barclay.argparser.Argument;

public class CollectF1R2CountsArgumentCollection {
    public static final String MIN_MEDIAN_MQ_LONG_NAME = "f1r2-median-mq";
    public static final String MIN_BASE_QUALITY_LONG_NAME = "f1r2-min-bq";
    public static final String MAX_DEPTH_LONG_NAME = "f1r2-max-depth";

    @Argument(fullName = MIN_MEDIAN_MQ_LONG_NAME, doc = "skip sites with median mapping quality below this value", optional = true)
    public int minMedianMapQual = 50;

    @Argument(fullName = MIN_BASE_QUALITY_LONG_NAME, doc = "exclude bases below this quality from pileup", optional = true)
    public int minBaseQuality = 20;

    @Argument(fullName = MAX_DEPTH_LONG_NAME, doc = "sites with depth higher than this value will be grouped", optional = true)
    public int maxDepth = F1R2FilterConstants.DEFAULT_MAX_DEPTH;
}
