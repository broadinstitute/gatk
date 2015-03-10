package org.broadinstitute.hellbender.tools.recalibration;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.covariates.ContextCovariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

public class ContextCovariateUnitTest {
    ContextCovariate covariate;
    RecalibrationArgumentCollection RAC;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new ContextCovariate();
        covariate.initialize(RAC);
    }

    @BeforeMethod
    public void initCache() {
        ReadCovariates.clearKeysCache();
    }

    @Test(enabled = true)
    public void testSimpleContexts() {
        SAMRecord read = ArtificialSAMUtils.createRandomRead(1000);
        SAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, RAC.LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);
        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);

        verifyCovariateArray(readCovariates.getMismatchesKeySet(), RAC.MISMATCHES_CONTEXT_SIZE, clippedRead, covariate);
        verifyCovariateArray(readCovariates.getInsertionsKeySet(), RAC.INDELS_CONTEXT_SIZE, clippedRead, covariate);
        verifyCovariateArray(readCovariates.getDeletionsKeySet(),  RAC.INDELS_CONTEXT_SIZE,  clippedRead, covariate);
    }

    public static void verifyCovariateArray(int[][] values, int contextSize, SAMRecord read, Covariate contextCovariate) {
        for (int i = 0; i < values.length; i++)
            Assert.assertEquals(contextCovariate.formatKey(values[i][0]), expectedContext(read, i, contextSize));

    }

    public static String expectedContext (SAMRecord read, int offset, int contextSize) {
        final String bases = stringFrom(read.getReadBases());
        String expectedContext = null;
        if (offset - contextSize + 1 >= 0) {
            String context = bases.substring(offset - contextSize + 1, offset + 1);
            if (!context.contains("N"))
                expectedContext = context;
        }
        return expectedContext;
    }

    private static String stringFrom(byte[] array) {
        String s = "";
        for (byte value : array)
            s += (char) value;
        return s;
    }

}
