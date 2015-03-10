package org.broadinstitute.hellbender.tools.recalibration;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.covariates.ReadGroupCovariate;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

public class ReadGroupCovariateUnitTest {
    ReadGroupCovariate covariate;
    RecalibrationArgumentCollection RAC;

    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        covariate = new ReadGroupCovariate();
        covariate.initialize(RAC);
    }

    @BeforeMethod
    public void initCache() {
        ReadCovariates.clearKeysCache();
    }

    @Test(enabled = true)
    public void testSingleRecord() {
        final String expected = "SAMPLE.1";
        SAMReadGroupRecord rg = new SAMReadGroupRecord("MY.ID");
        rg.setPlatformUnit(expected);
        runTest(rg, expected, covariate);
    }

    @Test(enabled = true)
    public void testMissingPlatformUnit() {
        final String expected = "MY.7";
        SAMReadGroupRecord rg = new SAMReadGroupRecord(expected);
        runTest(rg, expected, covariate);
    }

    @Test(enabled = true)
    public void testForceReadgroup() {
        final RecalibrationArgumentCollection forcedRAC = new RecalibrationArgumentCollection();
        forcedRAC.FORCE_READGROUP = "FOO";
        final ReadGroupCovariate forcedCovariate = new ReadGroupCovariate();
        forcedCovariate.initialize(forcedRAC);

        final SAMReadGroupRecord rg = new SAMReadGroupRecord("NOT_FOO");
        runTest(rg, "FOO", forcedCovariate);
    }

    private static void runTest(final SAMReadGroupRecord rg, final String expected, final ReadGroupCovariate covariate) {
        SAMRecord read = ArtificialSAMUtils.createRandomRead(10);
        ReadUtils.setReadGroup(read, rg);
        ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), 1);
        covariate.recordValues(read, readCovariates);
        verifyCovariateArray(readCovariates.getMismatchesKeySet(), expected, covariate);

    }

    private static void verifyCovariateArray(final int[][] values, final String expected, final ReadGroupCovariate covariate) {
        for (int[] value : values) {
            String actual = covariate.formatKey(value[0]);
            Assert.assertEquals(actual, expected);
        }
    }

}
