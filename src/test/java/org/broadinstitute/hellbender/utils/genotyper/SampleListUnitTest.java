package org.broadinstitute.hellbender.utils.genotyper;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Test {@link org.broadinstitute.gatk.utils.genotyper.AlleleListUtils}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleListUnitTest {

    @Test(dataProvider = "singleSampleListData")
    public void testAsList(final List<String> samples) {
         final SampleList sampleList = new IndexedSampleList(samples);
         final List<String> asList = sampleList.asListOfSamples();
         Assert.assertEquals(samples, asList);
    }

    @Test(dataProvider = "singleSampleListData")
    public void testAsSet(final List<String> samples) {
        final SampleList sampleList = new IndexedSampleList(samples);
        final Set<String> asSet = sampleList.asSetOfSamples();
        Assert.assertEquals(new LinkedHashSet<>(samples), asSet);
    }

    @Test
    public void testEmpty() {
        final SampleList sampleList = SampleList.emptySampleList();
        Assert.assertEquals(sampleList.numberOfSamples(), 0);
        Assert.assertTrue(sampleList.indexOfSample("bozo") < 0);

        final List<String> asList = sampleList.asListOfSamples();
        Assert.assertEquals(asList, Collections.emptyList());

        final Set<String> asSet = sampleList.asSetOfSamples();
        Assert.assertEquals(asSet, Collections.emptySet());
    }

        @Test
    public void testSingleton() {
        final String s = "fred";
        final SampleList sampleList = SampleList.singletonSampleList(s);
        final List<String> asList = sampleList.asListOfSamples();
        Assert.assertEquals(asList, Arrays.asList(s));
        Assert.assertEquals(sampleList.getSample(0), s);
        Assert.assertEquals(sampleList.indexOfSample(s), 0);
        Assert.assertNotEquals(sampleList.indexOfSample("bozo"), 0);

        Assert.assertTrue(asList.contains(s));
        Assert.assertTrue(!asList.contains("bozo"));

        final Set<String> asSet = sampleList.asSetOfSamples();
        Assert.assertEquals(asSet, new LinkedHashSet<>(Arrays.asList(s)));
        Assert.assertTrue(asSet.contains(s));
        Assert.assertTrue(! asSet.contains("bozo"));
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

  private List<String>[] sampleLists;

    @BeforeClass
    @SuppressWarnings({"unchecked", "rawtypes"})
    public void setUp() {
        sampleLists = (List<String>[])new List[SAMPLE_COUNT.length];
        int nextIndex = 0;
        for (int i = 0; i < SAMPLE_COUNT.length; i++) {
            final List<String> sampleList = new ArrayList<>(SAMPLE_COUNT[i]);
            sampleList.add("SAMPLE_" + i);
            sampleLists[nextIndex++] = sampleList;
        }
    }

    private static final int[] SAMPLE_COUNT = { 0, 1, 5, 10, 20};


    @DataProvider(name="singleSampleListData")
    public Object[][] singleSampleListData() {
        final Object[][] result = new Object[sampleLists.length][];
        for (int i = 0; i < sampleLists.length; i++)
            result[i] = new Object[] { sampleLists[i]};
        return result;
    }

    @DataProvider(name="twoSampleListData")
    public Object[][] twoAlleleListData() {
        final Object[][] result = new Object[sampleLists.length * sampleLists.length][];
        int index = 0;
        for (int i = 0; i < sampleLists.length; i++)
            for (int j = 0; j < sampleLists.length; j++)
                result[index++] = new Object[] { sampleLists[i], sampleLists[j]};
        return result;
    }







}
