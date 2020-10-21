package org.broadinstitute.hellbender.transformers;

import java.util.List;

import org.broadinstitute.hellbender.utils.read.GATKRead;

public final class BaseQualityStaticQuantizationTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;
    private byte[] staticQuantizedMapping;

    public BaseQualityStaticQuantizationTransformer(final List<Integer> staticQuantizationQuals, final boolean roundDown) {
        staticQuantizedMapping = BQSRReadTransformer.constructStaticQuantizedMapping(staticQuantizationQuals, roundDown);
    }

    /**
     * Performs static quantization of the quality scores
     *
     * @param originalRead the read to recalibrate
     */
    @Override
    public GATKRead apply(final GATKRead read) {
        final byte[] quals = read.getBaseQualities();
        final int readLength = quals.length;

        for (int offset = 0; offset < readLength; offset++) { // recalibrate all bases in the read
            quals[offset] = staticQuantizedMapping[quals[offset]];
        }
        read.setBaseQualities(quals);
        return read;
    }
}
