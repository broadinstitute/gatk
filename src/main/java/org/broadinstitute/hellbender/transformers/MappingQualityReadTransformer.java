package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A read transformer to modify the mapping quality of reads with MQ=255 to reads with MQ=60
 *
 *
 */
public class MappingQualityReadTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;

    private int fromQuality = 255;
    private int toQuality = 60;

    public MappingQualityReadTransformer(final int fromQuality, final int toQuality) {
        this.fromQuality = fromQuality;
        this.toQuality = toQuality;
    }

    @Override
    public GATKRead apply(final GATKRead read) {
        if (read.getMappingQuality() == fromQuality) {
            read.setMappingQuality(toQuality);
        }
        return read;
    }
}
