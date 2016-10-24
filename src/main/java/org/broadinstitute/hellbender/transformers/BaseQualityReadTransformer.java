package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.utils.read.GATKRead;

public final class BaseQualityReadTransformer implements ReadTransformer {

    private static final long serialVersionUID = 1L;

    private int BASE_QUALITY_THRESHOLD = 15;

    public BaseQualityReadTransformer( ) { }

    public BaseQualityReadTransformer(final int quality_threshold) {
        this.BASE_QUALITY_THRESHOLD = quality_threshold;
    }

    //Changes bases with qualities lower than BASE_QUALITY_THRESHOLD to 'N'
    @Override
    public GATKRead apply( final GATKRead read ) {
        if (read.getBaseQualityCount() == read.getLength()) {
            byte[] bases_new = read.getBases();
            for (int i = 0; i < read.getLength(); i++) {
                if (read.getBaseQuality(i) < BASE_QUALITY_THRESHOLD) {
                    bases_new[i] = 'N';
                }
            }
            read.setBases(bases_new);
        }
        return read;
    }
}
