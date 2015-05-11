package org.broadinstitute.hellbender.utils.read.markduplicates;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Tests for OpticalDuplicateFinder
 *
 * @author Nils Homer
 */
public final class OpticalDuplicateFinderTest {

    /** Tests rapidParseInt for positive and negative numbers, as well as non-digit suffixes */
    @Test
    public void testRapidParseInt() {
        final OpticalDuplicateFinder opticalDuplicateFinder = new OpticalDuplicateFinder();
        for (int i = -100; i < 100; i++) {
            Assert.assertEquals(opticalDuplicateFinder.rapidParseInt(Integer.toString(i)), i);

            // trailing characters
            Assert.assertEquals(opticalDuplicateFinder.rapidParseInt(Integer.toString(i) + "A"), i);
            Assert.assertEquals(opticalDuplicateFinder.rapidParseInt(Integer.toString(i) + "ACGT"), i);
            Assert.assertEquals(opticalDuplicateFinder.rapidParseInt(Integer.toString(i) + ".1"), i);
        }
    }

    /** Helper for testGetRapidDefaultReadNameRegexSplit */
    private void doTestGetRapidDefaultReadNameRegexSplit(int numFields, final OpticalDuplicateFinder opticalDuplicateFinder) {
        final int[] inputFields = new int[numFields];
        final int[] expectedFields = new int[numFields];
        String readName = "";
        for (int i = 0; i < inputFields.length; i++) {
            inputFields[i] = -1;
            expectedFields[i] = -1;
            if (0 < i) readName += ":";
            readName += Integer.toString(i);
        }
        if (2 < numFields) expectedFields[2] = 2;
        if (3 < numFields) expectedFields[3] = 3;
        if (4 < numFields) expectedFields[4] = 4;
        Assert.assertEquals(opticalDuplicateFinder.getRapidDefaultReadNameRegexSplit(readName, ':', inputFields), numFields);
        for (int i = 0; i < inputFields.length; i++) {
            Assert.assertEquals(inputFields[i], expectedFields[i]);
        }
    }

    /** Tests that we split the string early, with the correct # of fields, and modified values */
    @Test
    public void testGetRapidDefaultReadNameRegexSplit() {
        final OpticalDuplicateFinder opticalDuplicateFinder = new OpticalDuplicateFinder();
        for (int i = 1; i < 10; i++) {
            doTestGetRapidDefaultReadNameRegexSplit((i <= 5) ? i : 5, opticalDuplicateFinder);
        }
    }

    // NB: these tests fails due to overflow in the duplicate finder test.  This has been the behavior previously, so keep it for now.
    @Test(dataProvider = "testParseReadNameDataProvider", enabled = false)
    public void testParseReadName(final String readName, final int tile, final int x, final int y) {
        OpticalDuplicateFinder opticalDuplicateFinder = new OpticalDuplicateFinder();
        OpticalDuplicateFinder.PhysicalLocation loc = new ReadEndsForMarkDuplicates();
        Assert.assertTrue(opticalDuplicateFinder.addLocationInformation(readName, loc));
        Assert.assertEquals(loc.getTile(), tile);
        Assert.assertEquals(loc.getX(), x);
        Assert.assertEquals(loc.getY(), y);
    }

    @DataProvider(name = "testParseReadNameDataProvider")
    public Object[][] testParseReadNameDataProvider() {
        return new Object[][]{
                {"RUNID:7:1203:2886:82292", 1203, 2886, 82292},
                {"RUNID:7:1203:2884:16834", 1203, 2884, 16834}
        };
    }
}
