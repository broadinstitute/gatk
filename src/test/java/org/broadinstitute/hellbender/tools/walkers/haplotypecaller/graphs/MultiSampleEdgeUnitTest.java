package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class MultiSampleEdgeUnitTest extends BaseTest {

    private class MultiplicityTestProvider {
        final List<Integer> countsPerSample;
        final int numSamplesPruning;
        public MultiplicityTestProvider(final List<Integer> countsPerSample, final int numSamplesPruning) {
            this.countsPerSample = countsPerSample;
            this.numSamplesPruning = numSamplesPruning;
        }
    }

    @DataProvider(name = "MultiplicityData")
    public Object[][] makeMultiplicityData() {
        List<Object[]> tests = new ArrayList<>();

        final List<Integer> countsPerSample = Arrays.asList(0, 1, 2, 3, 4, 5);
        for ( final int numSamplesPruning : Arrays.asList(1, 2, 3) ) {
            for ( final int nSamples : Arrays.asList(1, 2, 3, 4, 5)) {
                for ( final List<Integer> perm : Utils.makePermutations(countsPerSample, nSamples, false) ) {
                    tests.add(new Object[]{new MultiplicityTestProvider(perm, numSamplesPruning)});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "MultiplicityData")
    public void testMultiplicity(final MultiplicityTestProvider cfg) {
        final MultiSampleEdge edge = new MultiSampleEdge(false, 0, cfg.numSamplesPruning);
        Assert.assertEquals(edge.getMultiplicity(), 0);
        Assert.assertEquals(edge.getPruningMultiplicity(), 0);

        final MultiSampleEdge copy = edge.copy();
        Assert.assertEquals(edge.getCurrentSingleSampleMultiplicity(), copy.getCurrentSingleSampleMultiplicity());
        Assert.assertEquals(edge.getDotLabel(), copy.getDotLabel());
        Assert.assertEquals(edge.getPruningMultiplicity(), copy.getPruningMultiplicity());
        Assert.assertEquals(edge.getMultiplicity(), copy.getMultiplicity());
        Assert.assertEquals(edge.getClass(), copy.getClass());

        int total = 0;
        for ( int i = 0; i < cfg.countsPerSample.size(); i++ ) {
            int countForSample = 0;
            for ( int count = 0; count < cfg.countsPerSample.get(i); count++ ) {
                edge.incMultiplicity(1);
                total++;
                countForSample++;
                Assert.assertEquals(edge.getMultiplicity(), total);
                Assert.assertEquals(edge.getCurrentSingleSampleMultiplicity(), countForSample);
            }
            edge.flushSingleSampleMultiplicity();
        }

        ArrayList<Integer> counts = new ArrayList<>(cfg.countsPerSample);
        counts.add(0);
        Collections.sort(counts);
        final int prune = counts.get(Math.max(counts.size() - cfg.numSamplesPruning, 0));
        Assert.assertEquals(edge.getMultiplicity(), total);
        Assert.assertEquals(edge.getPruningMultiplicity(), prune);
    }
}
