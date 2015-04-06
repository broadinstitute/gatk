package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Unit tests for {@link Sample}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleUnitTest extends BaseTest {

    /**
     * Different read-group counts used in testing.
     */
    private final int[] readGroupCount = {0,1,2,3,10,13,20,23,30,47,101};

    @Test
    public void testCreation() {
        new Sample("SM1");
    }

    @Test(dependsOnMethods = "testCreation")
    public void testReadGroupsOnRGLessSample() {
        final Sample sample = new Sample("SM1");
        Assert.assertNotNull(sample.readGroups());
        Assert.assertEquals(sample.readGroups().size(), 0);
    }

    @Test(dependsOnMethods = "testCreation")
    public void testSampleIdOnRGLessSample() {
        final Sample sample = new Sample("SM1");
        Assert.assertNotNull(sample.getId(),"SM1");
    }

    @Test
    public void testCreationSingleReadGroup() {
        new Sample("SM1", "RG1_1");
    }

    @Test(dependsOnMethods = "testCreationSingleReadGroup")
    public void testReadGroupsOnSingleRGSample() {
        final Sample sample = new Sample("SM1", "RG1_1");
        Assert.assertNotNull(sample.readGroups());
        Assert.assertEquals(sample.readGroups().size(),1);
        Assert.assertNotNull(sample.readGroups().iterator().next());
    }

    @Test(dependsOnMethods = "testCreationSingleReadGroup")
    public void testSampleIdOnSingleRGSample() {
        final Sample sample = new Sample("SM1", "RG1_1");
        Assert.assertEquals(sample.getId(),"SM1");
    }

    @Test(dataProvider = "multiReadGroupSampleData", dependsOnMethods = "testCreationSingleReadGroup")
    public void testSampleIdMultiReadGroup(final Collection<String> rgs) {
        final Sample sample = new Sample("SM1", rgs);
        Assert.assertEquals(sample.getId(),"SM1");
    }

    @Test(dataProvider = "multiReadGroupSampleData", dependsOnMethods = "testCreationSingleReadGroup")
    public void testReadGroupsMultiReadGroup(final Collection<String> rgs) {
        final Sample sample = new Sample("SM1", rgs);
        Assert.assertNotNull(sample.readGroups());
        Assert.assertEquals(sample.readGroups().size(),rgs.size());
        for (final String rg : rgs) {
            Assert.assertTrue(sample.readGroups().contains(rg), "missing read-group " + rg);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullIdCreation() {
        new Sample(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullReadGroupCreation() {
        new Sample("SM1",(String) null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    @SuppressWarnings("unchecked")
    public void testNullReadGroupCollectionCreation() {
        new Sample("SM1",(Collection)null);
    }

    @Test(expectedExceptions = {IllegalArgumentException.class})
    public void testNullReadGroupInCollection() {
        new Sample("SM1", Arrays.asList("RG0", "RG1", null, "RG3"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRepeatedReadGroupInCollection() {
        final Sample sample = new Sample("SM1", Arrays.asList("RG0", "RG1", "RG0"));
        Assert.assertEquals(sample.readGroups().size(),2);
    }

    @DataProvider(name="multiReadGroupSampleData")
    public Object[][] multiReadGroupSampleData() {
        final List<Object[]> result = new ArrayList<>(readGroupCount.length);
        for (final int rgCount : readGroupCount) {
            final List<String> rgs = new ArrayList<>(rgCount);
            for (int i = 0; i < rgCount; i++) {
                rgs.add("RG1_" + (i + 1));
            }
            result.add(new Object[] { rgs });
        }
        return result.toArray(new Object[result.size()][]);
    }
}
