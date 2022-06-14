package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;

public final class AltCovariate implements Covariate {
    private static final long serialVersionUID = 1L;

    public AltCovariate(final RecalibrationArgumentCollection RAC){
    }

    // Used to pick out the covariate's value from attributes of the read
    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final ReadCovariates values, final boolean recordIndelValues) {
        final int readLength = read.getLength();
        final byte[] readBases = read.getBasesNoCopy();
        for (int i = 0; i < readLength; i++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(readBases[i]);
            values.addCovariate(baseIndex, 0, 0, i);
        }
    }

    @Override
    public String formatKey(final int key){
        if ( key < 0 ) {
            return "N";
        } else {
            return String.format("%c", (char)BaseUtils.baseIndexToSimpleBase(key));
        }
    }

    @Override
    public int keyFromValue(final Object value) {
        if ( (value instanceof String)) {
            final byte base = ((String)value).getBytes()[0];
            if ( base == 'N' ) {
                return -1;
            } else {
                return BaseUtils.simpleBaseToBaseIndex(base);
            }

        } else {
            return (Integer)value;
        }
    }

    @Override
    public int maximumKeyValue() {
        return BaseUtils.Base.values().length;
    }
}