package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.walkers.genotyper.MinimalGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.jetbrains.annotations.NotNull;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

import static org.testng.Assert.*;

public class MultiPloidyGenotyperCacheUnitTest extends GATKBaseTest {


    @DataProvider
    public static Object[][] getGenotyperParams() {
        return new Object[][]{
                {"1:10-30", 2},
                {"1:70-85", 3},
                {"1:80-93", 3},
              //  {"1:89-92", 2},  TODO is this a bug?
                {"1:90-95", 1},
                {"1:85-120", 5},
                {"2:200", 2},
                {"2:1-500", 4}
        };
    }

    @Test(dataProvider = "getGenotyperParams")
    public void testGettingGenotypers(String interval, int expectedPloidy){
        final var detector = getOverlapDetector();
        final MultiPloidyGenotyperCache<MinimalGenotypingEngine> genotypingCache = new MultiPloidyGenotyperCache<>(this::getEngine, 2, detector);
        final MinimalGenotypingEngine genotypingEngine = genotypingCache.getGenotypingEngine(new SimpleInterval(interval));
        Assert.assertEquals(genotypingEngine.getGenotypeArgs().samplePloidy, expectedPloidy);
    }

    @NotNull
    private OverlapDetector<SimpleCount> getOverlapDetector() {
        final var counts = List.of(
                getCount("1:80", 3),
                getCount("1:90-95", 1),
                getCount("1:100-110", 5),
                getCount("2:85-90", 4));
        return OverlapDetector.create(counts);
    }

    @Test
    public void testReturnsTheSameGenotyper(){
        final var detector = getOverlapDetector();
        final MultiPloidyGenotyperCache<MinimalGenotypingEngine> genotypingCache = new MultiPloidyGenotyperCache<>(this::getEngine, 2, detector);
        final MinimalGenotypingEngine genotypingEngine = genotypingCache.getGenotypingEngine(new SimpleInterval("1:10-20"));
        final MinimalGenotypingEngine genotypingEngine2 = genotypingCache.getGenotypingEngine(new SimpleInterval("2:10-20"));
        Assert.assertSame(genotypingEngine, genotypingEngine2);
    }

    @Test
    public void testGetDefaultBeforeInit() {
        final var detector = getOverlapDetector();
        final MultiPloidyGenotyperCache<MinimalGenotypingEngine> genotypingCache = new MultiPloidyGenotyperCache<>(this::getEngine, 2, detector);
        Assert.assertNotNull(genotypingCache.getDefaultGenotypingEngine());
    }

    private MinimalGenotypingEngine getEngine(int ploidy){
        final var args = new StandardCallerArgumentCollection();
        args.genotypeArgs.samplePloidy = ploidy;
        final var samples = SampleList.singletonSampleList("sample1");
        return new MinimalGenotypingEngine(args, samples);
    }

    private static SimpleCount getCount(String location, int count){
        return new SimpleCount(new SimpleInterval(location), count);
    }

}