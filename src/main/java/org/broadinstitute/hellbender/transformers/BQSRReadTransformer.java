package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.BaseRecalibration;

import java.io.File;

public final class BQSRReadTransformer implements ReadTransformer{
    private final BaseRecalibration bqsr;

    public BQSRReadTransformer(File bqsrRecalFile, int quantizationLevels, boolean disableIndelQuals, final int preserveQLessThan, final boolean emitOriginalQuals, final double globalQScorePrior) {
        this.bqsr = new BaseRecalibration(bqsrRecalFile, quantizationLevels, disableIndelQuals, preserveQLessThan, emitOriginalQuals, globalQScorePrior);
    }

    @Override
    public SAMRecord apply(SAMRecord read) {
        bqsr.recalibrateRead(read);
        return read;
    }
}
