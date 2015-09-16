package org.broadinstitute.hellbender.utils.recalibration;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CycleCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public final class RecalibrationReportUnitTest extends BaseTest {
    @BeforeMethod
    public void init() {
        ReadCovariates.clearKeysCache();
    }

    private static RecalDatum createRandomRecalDatum(int maxObservations, int maxErrors) {
        final Random random = new Random();
        final int nObservations = random.nextInt(maxObservations);
        final int nErrors = Math.min(random.nextInt(maxErrors), nObservations);
        final int qual = random.nextInt(QualityUtils.MAX_SAM_QUAL_SCORE);
        return new RecalDatum((long)nObservations, (double)nErrors, (byte)qual);
    }

    @Test(expectedExceptions = UserException.class)
    public void testUnsupportedCovariates(){
        File file = new File(publicTestDir + "org/broadinstitute/hellbender/tools/" + "unsupported-covariates.table.gz");
        new RecalibrationReport(file);
    }

    @Test
    public void testOutput() {
        final int length = 100;

        List<Byte> quals = new ArrayList<>(QualityUtils.MAX_SAM_QUAL_SCORE + 1);
        List<Long> counts = new ArrayList<>(QualityUtils.MAX_SAM_QUAL_SCORE + 1);

        for (int i = 0;  i<= QualityUtils.MAX_SAM_QUAL_SCORE; i++) {
            quals.add((byte) i);
            counts.add(1L);
        }

        final QuantizationInfo quantizationInfo = new QuantizationInfo(quals, counts);
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        quantizationInfo.noQuantization();
        final String readGroupID = "id";
        final StandardCovariateList covariateList = new StandardCovariateList(RAC, Collections.singletonList(readGroupID));

        final SAMReadGroupRecord rg = new SAMReadGroupRecord(readGroupID);
        rg.setPlatform("illumina");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(rg);

        final GATKRead read = ArtificialReadUtils.createRandomRead(header, length, false);
        read.setReadGroup(rg.getReadGroupId());

        final byte [] readQuals = new byte[length];
        for (int i = 0; i < length; i++)
            readQuals[i] = 20;
        read.setBaseQualities(readQuals);

        final int expectedKeys = expectedNumberOfKeys(length, RAC.INDELS_CONTEXT_SIZE, RAC.MISMATCHES_CONTEXT_SIZE);
        int nKeys = 0;  // keep track of how many keys were produced
        final ReadCovariates rc = RecalUtils.computeCovariates(read, header, covariateList);

        final RecalibrationTables recalibrationTables = new RecalibrationTables(covariateList);
        final NestedIntegerArray<RecalDatum> rgTable = recalibrationTables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getQualityScoreTable();

        for (int offset = 0; offset < length; offset++) {

            for (EventType errorMode : EventType.values()) {

                final int[] covariates = rc.getKeySet(offset, errorMode);
                final int randomMax = errorMode == EventType.BASE_SUBSTITUTION ? 10000 : 100000;

                rgTable.put(createRandomRecalDatum(randomMax, 10), covariates[0], errorMode.ordinal());
                qualTable.put(createRandomRecalDatum(randomMax, 10), covariates[0], covariates[1], errorMode.ordinal());
                nKeys += 2;
                for (NestedIntegerArray<RecalDatum> covTable : recalibrationTables.getAdditionalTables()){
                    Covariate cov = recalibrationTables.getCovariateForTable(covTable);
                    final int covValue = covariates[covariateList.indexByClass(cov.getClass())];
                    if ( covValue >= 0 ) {
                        covTable.put(createRandomRecalDatum(randomMax, 10), covariates[0], covariates[1], covValue, errorMode.ordinal());
                        nKeys++;
                    }
                }
            }
        }
        Assert.assertEquals(nKeys, expectedKeys);
    }

    private static int expectedNumberOfKeys (int readLength, int indelContextSize, int mismatchesContextSize) {
        final int numCovariates = 4;
        final int numTables = 3;
        final int mismatchContextPadding = mismatchesContextSize - 1;
        final int indelContextPadding = 2 * (indelContextSize - 1);
        final int indelCyclePadding = 2 * (2 * CycleCovariate.CUSHION_FOR_INDELS);

        return (numCovariates * numTables * readLength) - mismatchContextPadding - indelContextPadding - indelCyclePadding;
    }

}
