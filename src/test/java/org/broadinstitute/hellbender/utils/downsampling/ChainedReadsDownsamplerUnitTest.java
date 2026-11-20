package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

public final class ChainedReadsDownsamplerUnitTest extends GATKBaseTest {

    private final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

    private List<GATKRead> stackedReads(final int lastStart, final int readsPerStart) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int start = 1; start <= lastStart; start++) {
            for (int i = 0; i < readsPerStart; i++) {
                reads.add(ArtificialReadUtils.createArtificialRead(header, "r" + start + "_" + i, 0, start, 100));
            }
        }
        return reads;
    }

    private static Map<Integer, Integer> depthByPosition(final List<GATKRead> reads) {
        final Map<Integer, Integer> depth = new TreeMap<>();
        for (final GATKRead read : reads) {
            for (int p = read.getStart(); p <= read.getEnd(); p++) {
                depth.merge(p, 1, Integer::sum);
            }
        }
        return depth;
    }

    @Test
    public void testBothStagesAreAppliedInOrder() {
        final List<GATKRead> reads = stackedReads(300, 10);
        final ChainedReadsDownsampler chained = new ChainedReadsDownsampler(
                new PositionalDownsampler(5, header, true),
                new MaxDepthDownsampler(50, 300, header, new Random(1)));
        Assert.assertTrue(chained.requiresCoordinateSortOrder());

        final List<GATKRead> out = new ArrayList<>();
        for (final GATKRead read : reads) {
            chained.submit(read);
            out.addAll(chained.consumeFinalizedItems());
        }
        chained.signalEndOfInput();
        out.addAll(chained.consumeFinalizedItems());

        // The positional stage leaves 5 per start (500x), which the depth stage then caps at 50.
        final Map<Integer, Integer> keptDepth = depthByPosition(out);
        Assert.assertTrue(keptDepth.get(150) >= 50 && keptDepth.get(150) < 200, "depth at 150 was " + keptDepth.get(150));
        Assert.assertEquals(chained.getNumberOfDiscardedItems() + out.size(), reads.size());
        Assert.assertEquals(new HashSet<>(out).size(), out.size(), "a read was emitted twice");
        for (int i = 1; i < out.size(); i++) {
            Assert.assertTrue(ReadCoordinateComparator.compareCoordinates(out.get(i - 1), out.get(i), header) <= 0);
        }
    }

    @Test
    public void testPendingAndFinalizedBookkeepingSpansBothStages() {
        final ChainedReadsDownsampler chained = new ChainedReadsDownsampler(
                new PositionalDownsampler(5, header, true),
                new MaxDepthDownsampler(50, 300, header, new Random(1)));
        Assert.assertFalse(chained.hasPendingItems());
        Assert.assertFalse(chained.hasFinalizedItems());
        Assert.assertEquals(chained.size(), 0);
        Assert.assertNull(chained.peekPending());
        Assert.assertNull(chained.peekFinalized());

        final GATKRead first = ArtificialReadUtils.createArtificialRead(header, "a", 0, 1, 100);
        chained.submit(first);
        // Held in the positional stage's reservoir until the position changes.
        Assert.assertTrue(chained.hasPendingItems());
        Assert.assertFalse(chained.hasFinalizedItems());
        Assert.assertEquals(chained.size(), 1);
        Assert.assertEquals(chained.peekPending(), first);

        chained.submit(ArtificialReadUtils.createArtificialRead(header, "b", 0, 2, 100));
        // "a" moved into the depth stage's window; still pending from the outside, and it is the pending read
        // nearest to being emitted, ahead of "b" which is still in the positional stage.
        Assert.assertTrue(chained.hasPendingItems());
        Assert.assertFalse(chained.hasFinalizedItems());
        Assert.assertEquals(chained.size(), 2);
        Assert.assertEquals(chained.peekPending(), first);

        chained.signalNoMoreReadsBefore(ArtificialReadUtils.createArtificialRead(header, "far", 0, 1000, 100));
        Assert.assertTrue(chained.hasFinalizedItems());
        Assert.assertEquals(chained.peekFinalized(), first);
        Assert.assertEquals(chained.consumeFinalizedItems().size(), 2);
        Assert.assertEquals(chained.size(), 0);
    }

    @Test
    public void testClearItemsAndResetStatsApplyToBothStages() {
        final ChainedReadsDownsampler chained = new ChainedReadsDownsampler(
                new PositionalDownsampler(1, header, true),
                new MaxDepthDownsampler(50, 300, header, new Random(1)));
        chained.submit(ArtificialReadUtils.createArtificialRead(header, "a", 0, 1, 100));
        chained.submit(ArtificialReadUtils.createArtificialRead(header, "b", 0, 1, 100));
        chained.submit(ArtificialReadUtils.createArtificialRead(header, "c", 0, 2, 100)); // finalizes position 1: one discarded
        Assert.assertEquals(chained.getNumberOfDiscardedItems(), 1);
        Assert.assertTrue(chained.hasPendingItems());

        chained.clearItems();
        chained.resetStats();

        Assert.assertFalse(chained.hasPendingItems());
        Assert.assertFalse(chained.hasFinalizedItems());
        Assert.assertEquals(chained.size(), 0);
        Assert.assertEquals(chained.getNumberOfDiscardedItems(), 0);
    }

    @DataProvider(name = "nullStages")
    public Object[][] nullStages() {
        final ReadsDownsampler stage = new MaxDepthDownsampler(50, 300, header, new Random(1));
        return new Object[][] {
                { null, stage },
                { stage, null },
        };
    }

    @Test(dataProvider = "nullStages", expectedExceptions = IllegalArgumentException.class)
    public void testNullStagesAreRejected(final ReadsDownsampler first, final ReadsDownsampler second) {
        new ChainedReadsDownsampler(first, second);
    }
}
