package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

public final class PassThroughDownsamplerUnitTest extends GATKBaseTest {
    @Test
    public void testClear() throws Exception {
        final PassThroughDownsampler ptd = new PassThroughDownsampler();
        for (int i = 0; i < 10; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead("100M");
            ptd.submit(read);
        }
        Assert.assertEquals(ptd.size(), 10);
        Assert.assertTrue(ptd.hasFinalizedItems());
        Assert.assertFalse(ptd.hasPendingItems());
        Assert.assertFalse(ptd.requiresCoordinateSortOrder());
        Assert.assertNotNull(ptd.peekFinalized());
        Assert.assertNull(ptd.peekPending());
        ptd.clearItems();
        Assert.assertEquals(ptd.size(), 0);
        Assert.assertFalse(ptd.hasFinalizedItems());
        Assert.assertFalse(ptd.hasPendingItems());
        Assert.assertFalse(ptd.requiresCoordinateSortOrder());
        Assert.assertNull(ptd.peekFinalized());
        Assert.assertNull(ptd.peekPending());

    }

    @Test
    public void testConsumeFinalizedItems() throws Exception {
        final PassThroughDownsampler ptd = new PassThroughDownsampler();
        for (int i = 0; i < 10; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead("100M");
            ptd.submit(read);
        }
        Assert.assertEquals(ptd.size(), 10);
        ptd.signalEndOfInput();
        ptd.signalNoMoreReadsBefore(ptd.peekFinalized());
        Assert.assertEquals(ptd.size(), 10);
        final List<GATKRead> finalizedItems = ptd.consumeFinalizedItems();
        Assert.assertEquals(finalizedItems.size(), 10);
        Assert.assertEquals(ptd.size(), 0);
    }


    @Test
    public void testSignalNoMoreReadsBefore() throws Exception {
        ReadsDownsampler rd = new PassThroughDownsampler();
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead("100M");
        final GATKRead r2 = ArtificialReadUtils.createArtificialRead("101M");
        rd.submit(r1);
        rd.signalNoMoreReadsBefore(r2);//no op
        rd.submit(r2);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoNullSignalNoMoreReadsBefore() throws Exception {
        ReadsDownsampler rd = new PassThroughDownsampler();
        rd.signalNoMoreReadsBefore(null);
    }
}