package org.broadinstitute.hellbender.utils.recalibration;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.covariates.ReadGroupCovariate;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class ReadGroupCovariateUnitTest {

    @BeforeMethod
    public void initCache() {
        ReadCovariates.clearKeysCache();
    }

    @Test
    public void testSingleRecord() {
        final String id = "MY.ID";
        final String expected = "SAMPLE.1";
        final ReadGroupCovariate covariate = new ReadGroupCovariate(new RecalibrationArgumentCollection(), Arrays.asList(expected));
        SAMReadGroupRecord rg = new SAMReadGroupRecord(id);
        rg.setPlatformUnit(expected);
        runTest(rg, expected, covariate);
    }

    @Test
    public void testMaxValue() {
        final String id = "MY.ID";
        final String expected = "SAMPLE.1";
        final ReadGroupCovariate covariate = new ReadGroupCovariate(new RecalibrationArgumentCollection(), Arrays.asList(expected));
        SAMReadGroupRecord rg = new SAMReadGroupRecord(id);
        rg.setPlatformUnit(expected);
        Assert.assertEquals(covariate.maximumKeyValue(), 0);//there's just 1 read group, so 0 is the max value
    }

    @Test
    public void testReadGroupNames() {
        final String id = "MY.ID";
        final String expected = "SAMPLE.1";
        final ReadGroupCovariate covariate = new ReadGroupCovariate(new RecalibrationArgumentCollection(), Arrays.asList(expected));
        final SAMFileHeader headerWithGroups = ArtificialReadUtils.createArtificialSamHeaderWithGroups(1, 0, 100, 2);
        final List<String> rgs = Arrays.asList("rg1", "rg2");
        Assert.assertEquals(ReadGroupCovariate.getReadGroupIDs(headerWithGroups), headerWithGroups.getReadGroups().stream().map(rg -> ReadGroupCovariate.getID(rg)).collect(Collectors.toList()));
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMissingKey() {
        final String id = "MY.ID";
        final String expected = "SAMPLE.1";
        final ReadGroupCovariate covariate = new ReadGroupCovariate(new RecalibrationArgumentCollection(), Arrays.asList(expected));
        final String s = covariate.formatKey(1);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMissingReadGroup() {
        final String id = "MY.ID";
        final String expected = "SAMPLE.1";
        final ReadGroupCovariate covariate = new ReadGroupCovariate(new RecalibrationArgumentCollection(), Arrays.asList(expected));
        final int key = covariate.keyFromValue("fred");
    }

    @Test
    public void testMissingPlatformUnit() {
        final String expected = "MY.7";
        final ReadGroupCovariate covariate = new ReadGroupCovariate(new RecalibrationArgumentCollection(), Arrays.asList(expected));
        SAMReadGroupRecord rg = new SAMReadGroupRecord(expected);
        runTest(rg, expected, covariate);
    }

    private static void runTest(final SAMReadGroupRecord rg, final String expected, final ReadGroupCovariate covariate) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(rg);

        GATKRead read = ArtificialReadUtils.createRandomRead(header, 10);
        read.setReadGroup(rg.getReadGroupId());

        ReadCovariates readCovariates = new ReadCovariates(read.getLength(), 1);
        covariate.recordValues(read, header, readCovariates);
        verifyCovariateArray(readCovariates.getMismatchesKeySet(), expected, covariate);

    }

    private static void verifyCovariateArray(final int[][] values, final String expected, final ReadGroupCovariate covariate) {
        for (int[] value : values) {
            String actual = covariate.formatKey(value[0]);
            Assert.assertEquals(actual, expected);
        }
    }

}
