package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleListUnitTester;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Tests {@link org.broadinstitute.gatk.tools.walkers.genotyper.HomogeneousPloidyModel}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HomogeneousPloidyModelUnitTest {
    private static final int[] PLOIDY = { 1, 2, 3, 7, 10};

    private static final int[] SAMPLE_COUNT = { 0, 1, 3, 4, 5, 6, 10, 101};


    @Test(dataProvider = "ploidyAndSampleListData")
    public void testPloidyAndSampleList(final int ploidy, final int sampleCount) {
        final List<String> sampleNames = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            sampleNames.add("SAMPLE_" + i);
        final IndexedSampleList sampleList = new IndexedSampleList(sampleNames);

        final HomogeneousPloidyModel ploidyModel = new HomogeneousPloidyModel(sampleList,ploidy);
        Assert.assertTrue(ploidyModel.isHomogeneous());
        Assert.assertEquals(ploidyModel.totalPloidy(), sampleCount * ploidy);

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
