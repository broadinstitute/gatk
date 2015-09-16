package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * The Reported Quality Score covariate.
 */
public final class QualityScoreCovariate implements Covariate {
    private static final long serialVersionUID = 1L;

    public QualityScoreCovariate(final RecalibrationArgumentCollection RAC){
        //nothing to initialize
    }

    @Override
    public void recordValues(final GATKRead read, final SAMFileHeader header, final ReadCovariates values) {
        final byte[] baseQualities = read.getBaseQualities();
        final byte[] baseInsertionQualities = ReadUtils.getBaseInsertionQualities(read);
        final byte[] baseDeletionQualities = ReadUtils.getBaseDeletionQualities(read);

        for (int i = 0; i < baseQualities.length; i++) {
            values.addCovariate(baseQualities[i], baseInsertionQualities[i], baseDeletionQualities[i], i);
        }
    }

    @Override
    public String formatKey(final int key) {
        return String.format("%d", key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return (value instanceof String) ? Byte.parseByte((String) value) : (int)(Byte) value;
    }

    @Override
    public int maximumKeyValue() {
        return QualityUtils.MAX_SAM_QUAL_SCORE;
    }
}