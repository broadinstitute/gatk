package org.broadinstitute.hellbender.utils.genotyper;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class SampleListUtilsUnitTest {

    private List<List<String>> sampleLists;

    @Test(dataProvider = "singleSampleListData")
    public void testAsList(final List<String> samples) {
         final SampleList sampleList = new IndexedSampleList(samples);
         final List<String> asList = sampleList.asListOfSamples();
         Assert.assertEquals(samples, asList);
    }

    @Test(dataProvider = "twoSampleListData", dependsOnMethods={"testAsList"})
    public void testEquals(final List<String> sample2, final List<String> samples2) {
        final SampleList sampleList1 = new IndexedSampleList(sample2);
        final SampleList sampleList2 = new IndexedSampleList(samples2);
        Assert.assertTrue(SampleList.equals(sampleList1, sampleList1));
        Assert.assertTrue(SampleList.equals(sampleList2, sampleList2));
        Assert.assertEquals(SampleList.equals(sampleList1, sampleList2),
                Arrays.equals(sampleList1.asListOfSamples().toArray(new String[sampleList1.numberOfSamples()]),
                              sampleList2.asListOfSamples().toArray(new String[sampleList2.numberOfSamples()]))
        );
        Assert.assertEquals(SampleList.equals(sampleList1, sampleList2),
                            SampleList.equals(sampleList2, sampleList1));
    }

    @BeforeClass
    public void setUp() {
        sampleLists = new ArrayList<>(SAMPLE_COUNT.length);
        for (int i = 0; i < SAMPLE_COUNT.length; i++) {
            sampleLists.add(i, null);
        }

        int nextIndex = 0;
        for (int i = 0; i < SAMPLE_COUNT.length; i++) {
            final List<String> sampleList = new ArrayList<>(SAMPLE_COUNT[i]);
            sampleList.add("SAMPLE_" + i);
            sampleLists.set(nextIndex++, sampleList);
        }
    }

    private static final int[] SAMPLE_COUNT = { 0, 1, 5, 10, 20};


    @DataProvider(name="singleSampleListData")
    public Object[][] singleSampleListData() {
        final Object[][] result = new Object[sampleLists.size()][];
        for (int i = 0; i < sampleLists.size(); i++)
            result[i] = new Object[] { sampleLists.get(i)};
        return result;
    }

    @DataProvider(name="twoSampleListData")
    public Object[][] twoAlleleListData() {
        final Object[][] result = new Object[sampleLists.size() * sampleLists.size()][];
        int index = 0;
        for (int i = 0; i < sampleLists.size(); i++)
            for (int j = 0; j < sampleLists.size(); j++)
                result[index++] = new Object[] { sampleLists.get(i), sampleLists.get(j)};
        return result;
    }







}
