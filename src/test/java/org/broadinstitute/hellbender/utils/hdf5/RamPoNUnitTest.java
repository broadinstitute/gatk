package org.broadinstitute.hellbender.utils.hdf5;


import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class RamPoNUnitTest extends BaseTest {

    private final static String TEST_DIR = "src/test/resources/org/broadinstitute/hellbender/tools/exome/";
    private final static String TEST_PCOV_FILE = "create-pon-control-full.pcov";

    @Test
    public void testCopiesReturnedOfMatrices() {
        final PoN filePoN = new HDF5PoN(new HDF5File(PoNTestUtils.createDummyHDF5FilePoN(new File(TEST_DIR + TEST_PCOV_FILE), 20), HDF5File.OpenMode.READ_ONLY));

        final PoN ramPoN = new RamPoN(filePoN);

        assertNormalizedCounts(filePoN, ramPoN);
        assertLogNormalizedCounts(filePoN, ramPoN);
        assertLogNormalizedPinvCounts(filePoN, ramPoN);
        assertReducedPanelCounts(filePoN, ramPoN);
        assertLogNormalizedPinvCounts(filePoN, ramPoN);
        assertReducedPanelPInverseCounts(filePoN, ramPoN);

        PoNTestUtils.assertEquivalentPoN(ramPoN, filePoN);
    }

    private void assertReducedPanelPInverseCounts(final PoN filePoN, final PoN ramPoN) {
        final RealMatrix ramNormalizedCounts = ramPoN.getReducedPanelPInverseCounts();
        final RealMatrix fileNormalizedCounts = filePoN.getReducedPanelPInverseCounts();

        Assert.assertEquals(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm(), 0, 1e-9);
        ramNormalizedCounts.setEntry(3, 4, 500000);
        Assert.assertEquals(ramNormalizedCounts.getEntry(3, 4), 500000, 1e-9);
        Assert.assertNotEquals(ramPoN.getReducedPanelPInverseCounts().getEntry(3, 4), 500000, 1e-9);

        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm() < 1e-9);

        final RealMatrix fileNormalizedCounts2 = filePoN.getReducedPanelPInverseCounts();
        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts2).getNorm() < 1e-9);
    }

    private void assertReducedPanelCounts(final PoN filePoN, final PoN ramPoN) {
        final RealMatrix ramNormalizedCounts = ramPoN.getReducedPanelCounts();
        final RealMatrix fileNormalizedCounts = filePoN.getReducedPanelCounts();

        Assert.assertEquals(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm(), 0, 1e-9);
        ramNormalizedCounts.setEntry(3, 4, 500000);
        Assert.assertEquals(ramNormalizedCounts.getEntry(3, 4), 500000, 1e-9);
        Assert.assertNotEquals(ramPoN.getReducedPanelCounts().getEntry(3, 4), 500000, 1e-9);

        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm() < 1e-9);

        final RealMatrix fileNormalizedCounts2 = filePoN.getReducedPanelCounts();
        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts2).getNorm() < 1e-9);
    }

    private void assertNormalizedCounts(final PoN filePoN, final PoN ramPoN) {
        final RealMatrix ramNormalizedCounts = ramPoN.getNormalizedCounts();
        final RealMatrix fileNormalizedCounts = filePoN.getNormalizedCounts();

        Assert.assertEquals(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm(), 0, 1e-9);
        ramNormalizedCounts.setEntry(3, 4, 500000);
        Assert.assertEquals(ramNormalizedCounts.getEntry(3, 4), 500000, 1e-9);
        Assert.assertNotEquals(ramPoN.getNormalizedCounts().getEntry(3, 4), 500000, 1e-9);

        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm() < 1e-9);

        final RealMatrix fileNormalizedCounts2 = filePoN.getNormalizedCounts();
        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts2).getNorm() < 1e-9);
    }

    private void assertLogNormalizedCounts(final PoN filePoN, final PoN ramPoN) {
        final RealMatrix ramNormalizedCounts = ramPoN.getLogNormalizedCounts();
        final RealMatrix fileNormalizedCounts = filePoN.getLogNormalizedCounts();

        Assert.assertEquals(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm(), 0, 1e-9);
        ramNormalizedCounts.setEntry(3, 4, 500000);
        Assert.assertEquals(ramNormalizedCounts.getEntry(3, 4), 500000, 1e-9);
        Assert.assertNotEquals(ramPoN.getLogNormalizedCounts().getEntry(3, 4), 500000, 1e-9);

        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm() < 1e-9);

        final RealMatrix fileNormalizedCounts2 = filePoN.getLogNormalizedCounts();
        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts2).getNorm() < 1e-9);
    }

    private void assertLogNormalizedPinvCounts(final PoN filePoN, final PoN ramPoN) {
        final RealMatrix ramNormalizedCounts = ramPoN.getLogNormalizedPInverseCounts();
        final RealMatrix fileNormalizedCounts = filePoN.getLogNormalizedPInverseCounts();

        Assert.assertEquals(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm(), 0, 1e-9);
        ramNormalizedCounts.setEntry(3, 4, 500000);
        Assert.assertEquals(ramNormalizedCounts.getEntry(3, 4), 500000, 1e-9);
        Assert.assertNotEquals(ramPoN.getLogNormalizedPInverseCounts().getEntry(3, 4), 500000, 1e-9);

        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts).getNorm() < 1e-9);

        final RealMatrix fileNormalizedCounts2 = filePoN.getLogNormalizedPInverseCounts();
        Assert.assertFalse(ramNormalizedCounts.subtract(fileNormalizedCounts2).getNorm() < 1e-9);
    }


}
