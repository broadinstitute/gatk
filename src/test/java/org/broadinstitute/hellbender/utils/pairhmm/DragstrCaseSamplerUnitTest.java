package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.IntervalPileup;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

public class DragstrCaseSamplerUnitTest {


    private static final String TEST_BAM = "/Users/valentin/Analysis/dragstr/NA12878.bam";
    private static final String TEST_REF = "/Users/valentin/Analysis/dragstr/download/DRAGstr/references/hg38_alt_aware.fa";

    @Test(dataProvider = "bugCasesData")
    public void testBugCases(final int chridx, final int pos, final String unit, final int repeats, final int expectedK, final int expectedN) {
        final ReferenceDataSource referenceDataSource = new ReferenceFileSource(new File(TEST_REF).toPath());
        final ReadsDataSource readsDataSource = new ReadsDataSource(new File(TEST_BAM).toPath());
        final DragstrModelEstimator estimator = new DragstrModelEstimator(new DragstrModelEstimatorArgumentCollection());
        final DragstrCasesSampler sampler = new DragstrCasesSampler(new DragstrCasesSamplerArgumentCollection(), referenceDataSource, readsDataSource);
        final DragstrLocus locus = DragstrLocus.make(chridx, pos, unit.getBytes(), repeats);
        final DragstrModelEstimator.RepeatCases cases = estimator.createPeriodCases(unit.length(), 20, 10).getRepeatCases(repeats);
        sampler.sample(cases, Collections.singletonList(locus));
        Assert.assertEquals(cases.size(), 1);
        Assert.assertEquals(cases.k[0], expectedK);
        Assert.assertEquals(cases.n[0], expectedN);
    }

    @DataProvider
    public Object[][] bugCasesData() {
        return new Object[][] {
                new Object[] {18, 28368088, "tttc", 18, 0, 1},
                new Object[] {1, 192311598, "tctt", 13, 1, 3},
                new Object[] {12, 22746119, "ta", 12, 1, 1},
                new Object[] {8, 73455283, "ccttc", 4, 0, 0}
        };
    }


}
