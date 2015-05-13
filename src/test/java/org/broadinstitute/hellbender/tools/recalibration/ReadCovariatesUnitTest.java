package org.broadinstitute.hellbender.tools.recalibration;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.tools.recalibration.covariates.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Random;

public final class ReadCovariatesUnitTest {

    @BeforeMethod
    public void init() {
        ReadCovariates.clearKeysCache();
    }

    @Test
    public void testCovariateGeneration() {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

        final String[] readGroups = {"RG1", "RG2", "RGbla"};
        ReadGroupCovariate rgCov = new ReadGroupCovariate(RAC, Arrays.asList(readGroups));
        QualityScoreCovariate qsCov = new QualityScoreCovariate(RAC);
        ContextCovariate coCov = new ContextCovariate(RAC);
        CycleCovariate cyCov = new CycleCovariate(RAC);

        StandardCovariateList covariates = new StandardCovariateList(RAC, Arrays.asList(readGroups));

        final int NUM_READS = 100;
        final Random rnd = Utils.getRandomGenerator();

        for (int idx = 0; idx < NUM_READS; idx++) {
            for (final String rgs : readGroups) {
                final int length = 10 + rnd.nextInt(100); // random read length, at least 10 bp long
                final SAMRecord read = ArtificialSAMUtils.createRandomRead(length, false);
                final SAMReadGroupRecord rg = new SAMReadGroupRecord(rgs);
                rg.setPlatform("illumina");
                ReadUtils.setReadGroup(read, rg);
                read.setReadNegativeStrandFlag(rnd.nextBoolean());
                final byte[] mQuals = read.getBaseQualities();
                final byte[] iQuals = ReadUtils.getBaseInsertionQualities(read);
                final byte[] dQuals = ReadUtils.getBaseDeletionQualities(read);
                ReadCovariates rc = RecalUtils.computeCovariates(read, covariates);

                // check that the length is correct
                Assert.assertEquals(rc.getMismatchesKeySet().length, length);
                Assert.assertEquals(rc.getInsertionsKeySet().length, length);
                Assert.assertEquals(rc.getDeletionsKeySet().length, length);

                for (int i = 0; i < length; i++) {
                    // check that read group is always the same
                    Assert.assertEquals(rgCov.formatKey(rc.getMismatchesKeySet(i)[0]), rgs);
                    Assert.assertEquals(rgCov.formatKey(rc.getInsertionsKeySet(i)[0]), rgs);
                    Assert.assertEquals(rgCov.formatKey(rc.getDeletionsKeySet(i)[0]), rgs);

                    // check quality score
                    Assert.assertEquals(qsCov.formatKey(rc.getMismatchesKeySet(i)[1]), String.valueOf(mQuals[i]));
                    Assert.assertEquals(qsCov.formatKey(rc.getInsertionsKeySet(i)[1]), String.valueOf(iQuals[i]));
                    Assert.assertEquals(qsCov.formatKey(rc.getDeletionsKeySet(i)[1]),  String.valueOf(dQuals[i]));

                    // check context
                    Assert.assertEquals(coCov.formatKey(rc.getMismatchesKeySet(i)[2]), ContextCovariateUnitTest.expectedContext(read, i, RAC.MISMATCHES_CONTEXT_SIZE, RAC.LOW_QUAL_TAIL), "read: " +idx  + " readGroup:" + rgs + " context mismatch key at position:" + i);
                    Assert.assertEquals(coCov.formatKey(rc.getInsertionsKeySet(i)[2]), ContextCovariateUnitTest.expectedContext(read, i, RAC.INDELS_CONTEXT_SIZE, RAC.LOW_QUAL_TAIL), "read: " +idx  + " readGroup:" + rgs + " context insertion key at position:" + i);
                    Assert.assertEquals(coCov.formatKey(rc.getDeletionsKeySet(i)[2]),  ContextCovariateUnitTest.expectedContext(read, i, RAC.INDELS_CONTEXT_SIZE, RAC.LOW_QUAL_TAIL), "read: " +idx  + " readGroup:" + rgs + " context deletion key at position:" + i);

                    // check cycle
                    final int expectedCycleMismatch = CycleCovariateUnitTest.expectedCycle(read, i, false, RAC.MAXIMUM_CYCLE_VALUE);
                    final int expectedCycleIndel = CycleCovariateUnitTest.expectedCycle(read, i, true, RAC.MAXIMUM_CYCLE_VALUE);
                    Assert.assertEquals(cyCov.formatKey(rc.getMismatchesKeySet(i)[3]), String.valueOf(expectedCycleMismatch), "read: " + idx + " cycle mismatch key at position:" + i);
                    Assert.assertEquals(cyCov.formatKey(rc.getInsertionsKeySet(i)[3]), String.valueOf(expectedCycleIndel),  "read: " +idx + " cycle insertion key at position:" + i);
                    Assert.assertEquals(cyCov.formatKey(rc.getDeletionsKeySet(i)[3]), String.valueOf(expectedCycleIndel),  "read: " +idx + " cycle deletion key at position:" + i);
                }

            }

        }

    }

}
