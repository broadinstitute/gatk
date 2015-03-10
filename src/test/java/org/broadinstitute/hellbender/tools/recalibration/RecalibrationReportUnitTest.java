package org.broadinstitute.hellbender.tools.recalibration;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.CycleCovariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class RecalibrationReportUnitTest {
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

    @Test
    public void testOutput() {
        final int length = 100;

        List<Byte> quals = new ArrayList<Byte>(QualityUtils.MAX_SAM_QUAL_SCORE + 1);
        List<Long> counts = new ArrayList<Long>(QualityUtils.MAX_SAM_QUAL_SCORE + 1);

        for (int i = 0;  i<= QualityUtils.MAX_SAM_QUAL_SCORE; i++) {
            quals.add((byte) i);
            counts.add(1L);
        }

        final QuantizationInfo quantizationInfo = new QuantizationInfo(quals, counts);
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        quantizationInfo.noQuantization();
        final StandardCovariateList requiredCovariates = new StandardCovariateList();

        final Covariate rgCovariate = requiredCovariates.getReadGroupCovariate();
        rgCovariate.initialize(RAC);

        final Covariate qsCovariate = requiredCovariates.getQualityScoreCovariate();
        qsCovariate.initialize(RAC);

        final Covariate cxCovariate = requiredCovariates.getContextCovariate();
        cxCovariate.initialize(RAC);

        final Covariate cyCovariate = requiredCovariates.getCycleCovariate();
        cyCovariate.initialize(RAC);

        final SAMReadGroupRecord rg = new SAMReadGroupRecord("id");
        rg.setPlatform("illumina");
        final SAMRecord read = ArtificialSAMUtils.createRandomRead(length, false);
        ReadUtils.setReadGroup(read, rg);
        final byte [] readQuals = new byte[length];
        for (int i = 0; i < length; i++)
            readQuals[i] = 20;
        read.setBaseQualities(readQuals);

        final int expectedKeys = expectedNumberOfKeys(length, RAC.INDELS_CONTEXT_SIZE, RAC.MISMATCHES_CONTEXT_SIZE);
        int nKeys = 0;                                                                                                  // keep track of how many keys were produced
        final ReadCovariates rc = RecalUtils.computeCovariates(read, requiredCovariates);

        final RecalibrationTables recalibrationTables = new RecalibrationTables(requiredCovariates);
        final NestedIntegerArray<RecalDatum> rgTable = recalibrationTables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> qualTable = recalibrationTables.getQualityScoreTable();

        for (int offset = 0; offset < length; offset++) {

            for (EventType errorMode : EventType.values()) {

                final int[] covariates = rc.getKeySet(offset, errorMode);
                final int randomMax = errorMode == EventType.BASE_SUBSTITUTION ? 10000 : 100000;

                rgTable.put(createRandomRecalDatum(randomMax, 10), covariates[0], errorMode.ordinal());
                qualTable.put(createRandomRecalDatum(randomMax, 10), covariates[0], covariates[1], errorMode.ordinal());
                nKeys += 2;
                for (int j = requiredCovariates.getOptionalCovariatesStartIndex(); j < requiredCovariates.size(); j++) {
                    final NestedIntegerArray<RecalDatum> covTable = recalibrationTables.getTable(j);
                    final int covValue = covariates[j];
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
