package org.broadinstitute.hellbender.utils.samples;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SampleUnitTest extends GATKBaseTest {

    @DataProvider(name="basicSamples")
    public Object[][] basicSamples() {
        return new Object[][] {
                { new Sample("1C", "fam1", "1M", "1F", Sex.UNKNOWN), "1C", "fam1", "1M", "1F", Sex.UNKNOWN, Affection.UNKNOWN },
                { new Sample("1F", "fam1", null, null, Sex.MALE), "1F", "fam1", null, null, Sex.MALE, Affection.UNKNOWN },
                { new Sample("1M", "fam1", null, null, Sex.FEMALE), "1M", "fam1", null, null, Sex.FEMALE, Affection.UNKNOWN },

                // Samples with Affection
                { new Sample("1C", "fam1", "1M", "1F", Sex.UNKNOWN, Affection.AFFECTED), "1C", "fam1", "1M", "1F", Sex.UNKNOWN, Affection.AFFECTED},
                { new Sample("1F", null, null, null, Sex.MALE, Affection.UNAFFECTED), "1F", null, null, null, Sex.MALE, Affection.UNAFFECTED },
                { new Sample("1M", null, null, null, Sex.FEMALE, Affection.OTHER), "1M", null, null, null, Sex.FEMALE, Affection.OTHER }
        };
    }

    /**
     * Basic getters
     */
    @Test(dataProvider="basicSamples")
    public void basicSampleTest(Sample sample, String id, String famID, String paternalID, String maternalID, Sex gender, Affection affection) {
        Assert.assertTrue(id.equals(sample.getID()));
        Assert.assertTrue(famID == null || famID.equals(sample.getFamilyID()));
        Assert.assertTrue(maternalID == null || maternalID.equals(sample.getMaternalID()));
        Assert.assertTrue(paternalID == null || paternalID.equals(sample.getPaternalID()));
        Assert.assertEquals(gender, sample.getSex());
        Assert.assertEquals(affection, sample.getAffection());
    }

    @Test(dataProvider="basicSamples")
    public void testMergeSamples(Sample sample, String id, String famID, String paternalID, String maternalID, Sex gender, Affection affection) {

        Sample newSample = new Sample("newSample", null, null, null, Sex.UNKNOWN, Affection.UNKNOWN);
        Sample mergedSample1 = newSample.mergeSamples(sample);
        Assert.assertTrue(mergedSample1.getID().equals("newSample"));

        if (famID == null) {
            Assert.assertEquals(null, mergedSample1.getFamilyID());
        }
        else {
            Assert.assertTrue(famID.equals(mergedSample1.getFamilyID()));
        }

        if (maternalID == null) {
            Assert.assertEquals(null, mergedSample1.getMaternalID());
        }
        else {
            Assert.assertTrue(maternalID.equals(mergedSample1.getMaternalID()));
        }

        if (paternalID == null) {
            Assert.assertEquals(null, mergedSample1.getPaternalID());
        }
        else {
            Assert.assertTrue(paternalID.equals(mergedSample1.getPaternalID()));
        }

        Assert.assertEquals(mergedSample1.getSex(), gender);
        Assert.assertEquals(mergedSample1.getAffection(), affection);
    }

    @DataProvider(name="sortSamplesNullFields")
    private Object[][] sortSamplesNullFields() {
        return new Object[][] {
                { new Sample("A", null, null, null, Sex.MALE, Affection.UNKNOWN), 0 },
                { new Sample("A", "fam1", null, null, Sex.MALE, Affection.UNKNOWN), -1 },
                { new Sample("A", null, "1M", null, Sex.MALE, Affection.UNKNOWN), -1 },
                { new Sample("A", null, null, "1F", Sex.MALE, Affection.UNKNOWN), -1 },
                { new Sample("A", null, null, null, Sex.FEMALE, Affection.UNKNOWN), -1 },
                { new Sample("A", null, null, null, Sex.MALE, Affection.AFFECTED), -1 }
        };
    }

    @Test(dataProvider="sortSamplesNullFields")
    private void testCompareNullFields(Sample target, int expected) {
        Sample source = new Sample("A", null, null, null, Sex.MALE, Affection.UNKNOWN);
        Assert.assertEquals(source.compareTo(target), expected);
    }

    @DataProvider(name="sortSamplesFields")
    private Object[][] sortSamplesFields() {
        return new Object[][] {
                { new Sample("A", "fam1", "1M", "1F", Sex.MALE, Affection.AFFECTED), 0 },
                { new Sample("Z", "fam1", "1M", "1F", Sex.MALE, Affection.AFFECTED), -25 },
                { new Sample("A", "fam2", "1M", "1F", Sex.MALE, Affection.AFFECTED), -1 },
                { new Sample("A", "fam1", "2M", "1F", Sex.MALE, Affection.AFFECTED), -1 },
                { new Sample("A", "fam1", "1M", "2F", Sex.MALE, Affection.AFFECTED), -1 },
                { new Sample("A", "fam1", "1M", "1F", Sex.FEMALE, Affection.AFFECTED), -1 },
                { new Sample("A", "fam1", "1M", "1F", Sex.MALE, Affection.UNAFFECTED), -1 },
        };
    }

    @Test(dataProvider="sortSamplesFields")
    private void testCompareFields(Sample target, int expected) {
        Sample source = new Sample("A", "fam1", "1M", "1F", Sex.MALE, Affection.AFFECTED);
        Assert.assertEquals(source.compareTo(target), expected);
    }

}
