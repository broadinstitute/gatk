package org.broadinstitute.hellbender.tools.recalibration.covariates;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;

/**
 * The Reported Quality Score covariate.
 */

public class QualityScoreCovariate implements Covariate {

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {}

    @Override
    public void recordValues(final SAMRecord read, final ReadCovariates values) {
        final byte[] baseQualities = read.getBaseQualities();
        final byte[] baseInsertionQualities = ReadUtils.getBaseInsertionQualities(read);
        final byte[] baseDeletionQualities = ReadUtils.getBaseDeletionQualities(read);

        for (int i = 0; i < baseQualities.length; i++) {
            values.addCovariate((int)baseQualities[i], (int)baseInsertionQualities[i], (int)baseDeletionQualities[i], i);
        }
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return Byte.parseByte(str);
    }

    @Override
    public String formatKey(final int key) {
        return String.format("%d", key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return (value instanceof String) ? (int) Byte.parseByte((String) value) : (int)(Byte) value;
    }

    @Override
    public int maximumKeyValue() {
        return QualityUtils.MAX_SAM_QUAL_SCORE;
    }
}