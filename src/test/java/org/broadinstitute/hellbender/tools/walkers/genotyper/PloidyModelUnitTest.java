package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleListUnitTester;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Tests {@link org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class PloidyModelUnitTest {
    private static final int[][] PLOIDIES =
            {{1, 1, 1, 1},
             {2, 2, 2, 2},
             {2, 4, 2, 4, 5, 6},
             {1, 2, 3, 7, 10}
            };


    @Test(dataProvider = "ploidiesAndSampleListData")
    public void testPloidyAndSampleList(final int[] ploidies) {
        final int sampleCount = ploidies.length;
        final List<String> sampleNames = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            sampleNames.add("SAMPLE_" + i);
        final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

        final PloidyModel ploidyModel = new PloidyModel(sampleList,ploidies);
        final boolean expectedHom = allSame(ploidies);


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

        new PloidyModel(sampleList, new int[]{1,2,3});//count mismatch

    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadPloidy() throws Exception {
        final int sampleCount = 2;
        final List<String> sampleNames = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            sampleNames.add("SAMPLE_" + i);
        final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

        new PloidyModel(sampleList, new int[]{1,-2});//bad ploidy

    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadSampleIndex() throws Exception {
        final int sampleCount = 2;
        final List<String> sampleNames = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            sampleNames.add("SAMPLE_" + i);
        final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

        final PloidyModel model = new PloidyModel(sampleList, new int[]{1, 2});
        model.samplePloidy(-3);
    }

    @DataProvider(name="ploidiesAndSampleListData")
    public Object[][] ploidiesAndSampleListData() {
        final Object[][] result = new Object[PLOIDIES.length][];
        int index = 0;
        for (int i = 0; i < PLOIDIES.length; i++) {
            result[index++] = new Object[]{PLOIDIES[i]};
        }
        return result;
    }

















        private static final int[] PLOIDY = { 1, 2, 3, 7, 10};

        private static final int[] SAMPLE_COUNT = { 0, 1, 3, 4, 5, 6, 10, 101};


        @Test(dataProvider = "ploidyAndSampleListData")
        public void testPloidyAndSampleList(final int ploidy, final int sampleCount) {
            final List<String> sampleNames = new ArrayList<>(sampleCount);
            for (int i = 0; i < sampleCount; i++)
                sampleNames.add("SAMPLE_" + i);
            final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

            final PloidyModel ploidyModel = new PloidyModel(sampleList,ploidy);



            for (int i = 0; i < sampleCount; i++)
                Assert.assertEquals(ploidyModel.samplePloidy(i), ploidy);

            SampleListUnitTester.assertSampleList(ploidyModel, sampleNames);
        }

        @DataProvider(name="ploidyAndSampleListData")
        public Object[][] ploidyAndSampleListData() {
            final Object[][] result = new Object[PLOIDY.length * SAMPLE_COUNT.length][];
            int index = 0;
            for (int i = 0; i < PLOIDY.length; i++)
                for (int j = 0; j < SAMPLE_COUNT.length; j++ )
                    result[index++] = new Object[] { PLOIDY[i], SAMPLE_COUNT[j]};
            return result;
        }




}




