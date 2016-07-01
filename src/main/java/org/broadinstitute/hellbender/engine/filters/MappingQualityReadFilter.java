package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * Keep reads with mapping qualities above a specified threshold.
 */
public final class MappingQualityReadFilter implements ReadFilter, CommandLineFilter, Serializable {
    private static final long serialVersionUID = 1L;

    private final static String mappingQualityArgName = "mappingQuality";

    @Argument(fullName=mappingQualityArgName, shortName="mq", optional=true)
    public int minMappingQualityScore = 10;

    public MappingQualityReadFilter( final int minMappingQualityScore ) {
        this.minMappingQualityScore = minMappingQualityScore;
    }

    @Override
    public boolean test( final GATKRead read ) {
        return read.getMappingQuality() >= minMappingQualityScore;
    }

    @Override
    public String validate() {
        String message = null;
        if (minMappingQualityScore <= 0 || minMappingQualityScore > 255) {
            message = "requires a value in the range [0, 255] for \"" + mappingQualityArgName + "\"";
        }

        return message;
    }

}