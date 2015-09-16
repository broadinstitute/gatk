package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleListUnitTester;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Tests {@link org.broadinstitute.gatk.tools.walkers.genotyper.HeterogeneousPloidyModel}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HeterogeneousPloidyModelUnitTest {
    private static final int[][] PLOIDIES =
            {{1, 1, 1, 1},
             {2, 2, 2, 2},
             {2, 4, 2, 4, 5, 6},
             {1, 2, 3, 7, 10}
            };


    @Test(dataProvider = "ploidyAndSampleListData")
    public void testPloidyAndSampleList(final int[] ploidies) {
        final int sampleCount = ploidies.length;
        final List<String> sampleNames = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            sampleNames.add("SAMPLE_" + i);
        final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

        final PloidyModel ploidyModel = new HeterogeneousPloidyModel(sampleList,ploidies);
        final boolean expectedHom = allSame(ploidies);
        Assert.assertEquals(ploidyModel.isHomogeneous(), expectedHom);
        Assert.assertEquals(ploidyModel.totalPloidy(), MathUtils.sum(ploidies));

        for (int i = 0; i < sampleCount; i++) {
            Assert.assertEquals(ploidyModel.samplePloidy(i), ploidies[i]);
        }

        SampleListUnitTester.assertSampleList(ploidyModel, sampleNames);
    }

    private boolean allSame(int[] ploidies) {
        return IntStream.of(ploidies).allMatch(p -> p == ploidies[0]);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadLength() throws Exception {
        final int sampleCount = 2;
        final List<String> sampleNames = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            sampleNames.add("SAMPLE_" + i);
        final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

        new HeterogeneousPloidyModel(sampleList, new int[]{1,2,3});//count mismatch

    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadPloidy() throws Exception {
        final int sampleCount = 2;
        final List<String> sampleNames = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            sampleNames.add("SAMPLE_" + i);
        final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

        new HeterogeneousPloidyModel(sampleList, new int[]{1,-2});//bad ploidy

    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadSampleIndex() throws Exception {
        final int sampleCount = 2;
        final List<String> sampleNames = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            sampleNames.add("SAMPLE_" + i);
        final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

        final PloidyModel model = new HeterogeneousPloidyModel(sampleList, new int[]{1, 2});
        model.samplePloidy(-3);
    }

    @DataProvider(name="ploidyAndSampleListData")
    public Object[][] ploidyAndSampleListData() {
        final Object[][] result = new Object[PLOIDIES.length][];
        int index = 0;
        for (int i = 0; i < PLOIDIES.length; i++) {
            result[index++] = new Object[]{PLOIDIES[i]};
        }
        return result;
    }
}
