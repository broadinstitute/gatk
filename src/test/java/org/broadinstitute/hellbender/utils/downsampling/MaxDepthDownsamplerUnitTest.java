package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

public final class MaxDepthDownsamplerUnitTest extends GATKBaseTest {

    private static final int READ_LENGTH = 100;

    private static SAMFileHeader headerWithSample(final String sample) {
        final SAMReadGroupRecord readGroup = new SAMReadGroupRecord("rg-" + sample);
        readGroup.setSample(sample);
        return ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(readGroup);
    }

    private static GATKRead read(final SAMFileHeader header, final String name, final int start) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, start, READ_LENGTH);
        if (!header.getReadGroups().isEmpty()) {
            read.setReadGroup(header.getReadGroups().get(0).getId());
        }
        return read;
    }

    /** Coordinate-sorted reads with {@code depthPerStart} reads at every start in [firstStart, lastStart]. */
    private static List<GATKRead> pileOfReads(final SAMFileHeader header, final int firstStart, final int lastStart, final int depthPerStart) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int start = firstStart; start <= lastStart; start++) {
            for (int i = 0; i < depthPerStart; i++) {
                reads.add(read(header, "r" + start + "_" + i, start));
            }
        }
        return reads;
    }

    private static List<GATKRead> runDownsampler(final ReadsDownsampler downsampler, final List<GATKRead> reads) {
        final List<GATKRead> out = new ArrayList<>();
        for (final GATKRead read : reads) {
            downsampler.submit(read);
            out.addAll(downsampler.consumeFinalizedItems());
        }
        downsampler.signalEndOfInput();
        out.addAll(downsampler.consumeFinalizedItems());
        return out;
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

    private static Set<String> names(final List<GATKRead> reads) {
        return reads.stream().map(GATKRead::getName).collect(Collectors.toSet());
    }

    /** Asserts the downsampler's guarantee: every position keeps at least min(original depth, cap) reads. */
    private static void assertDepthFloorHolds(final List<GATKRead> original, final List<GATKRead> kept, final int cap) {
        final Map<Integer, Integer> originalDepth = depthByPosition(original);
        final Map<Integer, Integer> keptDepth = depthByPosition(kept);
        for (final Map.Entry<Integer, Integer> entry : originalDepth.entrySet()) {
            final int floor = Math.min(entry.getValue(), cap);
            final int actual = keptDepth.getOrDefault(entry.getKey(), 0);
            Assert.assertTrue(actual >= floor, "position " + entry.getKey() + " kept " + actual + " reads, below the floor of " + floor);
        }
    }

    private static void assertCoordinateSorted(final List<GATKRead> reads, final SAMFileHeader header) {
        for (int i = 1; i < reads.size(); i++) {
            Assert.assertTrue(ReadCoordinateComparator.compareCoordinates(reads.get(i - 1), reads.get(i), header) <= 0,
                    "output not coordinate sorted at index " + i);
        }
    }

    @Test
    public void testReadsBelowTheCapPassThroughUnchanged() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = pileOfReads(header, 1, 500, 1); // depth peaks at READ_LENGTH = 100
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(100, 300, header, new Random(1));
        Assert.assertTrue(downsampler.requiresCoordinateSortOrder());

        final List<GATKRead> out = runDownsampler(downsampler, reads);

        Assert.assertEquals(out, reads);
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);
    }

    @DataProvider(name = "capsAndWindows")
    public Object[][] capsAndWindows() {
        return new Object[][] {
                // cap, window width, reads per start, whether the window is wide enough for the kept depth to stay near the cap
                { 150, 300, 20, true },
                { 150, 1000, 20, true },
                { 150, 300, 2, true },    // original depth (200x) only modestly above the cap
                { 1, 300, 5, true },
                { 50, 50, 20, false },    // window narrower than a read: the floor holds but the cap overshoots
        };
    }

    @Test(dataProvider = "capsAndWindows")
    public void testDepthOverTheCapIsReducedButNeverBelowTheFloor(final int cap, final int windowSize, final int readsPerStart, final boolean windowIsWiderThanReads) {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = pileOfReads(header, 1, 1000, readsPerStart);
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(cap, windowSize, header, new Random(7));

        final List<GATKRead> out = runDownsampler(downsampler, reads);

        assertDepthFloorHolds(reads, out, cap);
        final int interiorDepth = depthByPosition(out).get(500);
        if (windowIsWiderThanReads) {
            // Reads are judged before later reads covering their right ends arrive, so kept depth overshoots the
            // cap by about one read's worth at window edges; well inside a window it stays close to the cap.
            Assert.assertTrue(interiorDepth < 2 * cap + READ_LENGTH, "kept depth in the interior was " + interiorDepth + " for a cap of " + cap);
        }
        Assert.assertTrue(out.size() < reads.size(), "nothing was discarded");
        Assert.assertEquals(out.size() + downsampler.getNumberOfDiscardedItems(), reads.size());
        assertCoordinateSorted(out, header);
    }

    @Test
    public void testPositionsUnderTheCapKeepEveryReadEvenNextToDeepPositions() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = new ArrayList<>();
        reads.addAll(pileOfReads(header, 1, 200, 1));      // shallow flank
        reads.addAll(pileOfReads(header, 201, 400, 10));   // deep block
        reads.addAll(pileOfReads(header, 401, 600, 1));    // shallow flank
        reads.sort(new ReadCoordinateComparator(header));
        final int cap = 120;
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(cap, 300, header, new Random(3));

        final List<GATKRead> out = runDownsampler(downsampler, reads);

        final Set<String> keptNames = names(out);
        final Map<Integer, Integer> originalDepth = depthByPosition(reads);
        for (final GATKRead read : reads) {
            boolean touchesShallowPosition = false;
            for (int p = read.getStart(); p <= read.getEnd() && !touchesShallowPosition; p++) {
                touchesShallowPosition = originalDepth.get(p) < cap;
            }
            if (touchesShallowPosition) {
                Assert.assertTrue(keptNames.contains(read.getName()), "read " + read.getName() + " touching a position under the cap was discarded");
            }
        }
        Assert.assertTrue(downsampler.getNumberOfDiscardedItems() > 0, "the deep block should have been capped");
        assertDepthFloorHolds(reads, out, cap);
    }

    @Test
    public void testDeletionsCountTowardsDepthAcrossTheirWholeSpan() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = new ArrayList<>();
        // 10 reads per start, each spanning 150 reference bases with a 50 base deletion in the middle.
        for (int start = 1; start <= 300; start++) {
            for (int i = 0; i < 10; i++) {
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "d" + start + "_" + i, 0, start,
                        new byte[READ_LENGTH], new byte[READ_LENGTH], "50M50D50M");
                read.setReadGroup(header.getReadGroups().get(0).getId());
                reads.add(read);
            }
        }
        final int cap = 100;

        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(cap, 300, header, new Random(11)), reads);

        assertDepthFloorHolds(reads, out, cap);
        Assert.assertTrue(out.size() < reads.size());
    }

    @Test
    public void testFloorHoldsAcrossWindowBoundaries() {
        final SAMFileHeader header = headerWithSample("s1");
        // Uniformly deep reads with windows narrower than a read, so every read spans at least one boundary and
        // the kept reads carried over from earlier windows are what keep later positions at the floor.
        final List<GATKRead> reads = pileOfReads(header, 1, 400, 10);
        final int cap = 50;

        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(cap, 60, header, new Random(5)), reads);

        assertDepthFloorHolds(reads, out, cap);
        final Map<Integer, Integer> keptDepth = depthByPosition(out);
        for (final int boundary : new int[]{61, 121, 181, 241, 301, 361}) {
            Assert.assertTrue(keptDepth.get(boundary) >= cap, "position " + boundary + " just past a window boundary kept only " + keptDepth.get(boundary));
        }
        Assert.assertTrue(out.size() < reads.size());
    }

    @Test
    public void testStateResetsAcrossContigs() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final List<GATKRead> reads = new ArrayList<>();
        for (final int contigIndex : new int[]{0, 1}) {
            for (int start = 1; start <= 200; start++) {
                for (int i = 0; i < 5; i++) {
                    reads.add(ArtificialReadUtils.createArtificialRead(header, "c" + contigIndex + "_" + start + "_" + i, contigIndex, start, READ_LENGTH));
                }
            }
        }
        final int cap = 100;

        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(cap, 300, header, new Random(9)), reads);

        for (final int contigIndex : new int[]{0, 1}) {
            final String contig = header.getSequence(contigIndex).getSequenceName();
            final List<GATKRead> original = reads.stream().filter(r -> r.getContig().equals(contig)).collect(Collectors.toList());
            final List<GATKRead> kept = out.stream().filter(r -> r.getContig().equals(contig)).collect(Collectors.toList());
            assertDepthFloorHolds(original, kept, cap);
            Assert.assertTrue(kept.size() < original.size(), "contig " + contig + " was not capped");
            // The first reads on each contig start on an empty depth array, so the first cap reads at position 1 are kept.
            Assert.assertTrue(depthByPosition(kept).get(1) >= Math.min(5, cap));
        }
        assertCoordinateSorted(out, header);
    }

    @Test
    public void testDifferentSeedsKeepDifferentSubsetsAndTheSameSeedIsReproducible() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = pileOfReads(header, 1, 600, 10);

        final Set<String> first = names(runDownsampler(new MaxDepthDownsampler(100, 300, header, new Random(1)), reads));
        final Set<String> second = names(runDownsampler(new MaxDepthDownsampler(100, 300, header, new Random(2)), reads));
        final Set<String> firstAgain = names(runDownsampler(new MaxDepthDownsampler(100, 300, header, new Random(1)), reads));

        Assert.assertNotEquals(first, second, "different seeds should keep different subsets");
        Assert.assertEquals(first, firstAgain, "the same seed must keep the same subset");
    }

    @Test
    public void testArrivalOrderModeIsDeterministicAndKeepsTheEarliestReads() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = pileOfReads(header, 1, 300, 10);

        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(50, 300, header, null), reads);
        final List<GATKRead> outAgain = runDownsampler(new MaxDepthDownsampler(50, 300, header, null), reads);

        Assert.assertEquals(out, outAgain);
        // In arrival order the first 50 reads (5 starts x 10 reads) are always kept.
        for (int i = 0; i < 50; i++) {
            Assert.assertEquals(out.get(i).getName(), reads.get(i).getName());
        }
    }

    @Test
    public void testDepthIsCappedPerSample() {
        final SAMReadGroupRecord readGroup1 = new SAMReadGroupRecord("rg1");
        readGroup1.setSample("s1");
        final SAMReadGroupRecord readGroup2 = new SAMReadGroupRecord("rg2");
        readGroup2.setSample("s2");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.addReadGroup(readGroup1);
        header.addReadGroup(readGroup2);
        final List<GATKRead> reads = new ArrayList<>();
        for (int start = 1; start <= 300; start++) {
            for (int i = 0; i < 4; i++) {
                final GATKRead deepSampleRead = ArtificialReadUtils.createArtificialRead(header, "s1_" + start + "_" + i, 0, start, READ_LENGTH);
                deepSampleRead.setReadGroup("rg1");
                reads.add(deepSampleRead);
            }
            final GATKRead shallowSampleRead = ArtificialReadUtils.createArtificialRead(header, "s2_" + start, 0, start, READ_LENGTH);
            shallowSampleRead.setReadGroup("rg2");
            reads.add(shallowSampleRead);
        }

        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(150, 300, header, new Random(5)), reads);

        final List<GATKRead> shallowSampleKept = out.stream().filter(r -> r.getReadGroup().equals("rg2")).collect(Collectors.toList());
        // Sample 2 peaks at 100x, under the cap, so all of its reads survive regardless of sample 1's depth.
        Assert.assertEquals(shallowSampleKept.size(), 300);
        Assert.assertTrue(out.size() < reads.size(), "sample 1 (400x) should have been capped");
    }

    @Test
    public void testReadsWithoutAnAssignedPositionPassThroughAfterBufferedReads() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = new ArrayList<>(pileOfReads(header, 1, 50, 1));
        final GATKRead unmapped = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30});
        reads.add(unmapped);

        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(1000, 300, header, new Random(1)), reads);

        Assert.assertEquals(out, reads, "the unplaced read must come out after the positioned reads it followed");
    }

    @Test
    public void testDownsamplerCanBeReusedAfterEndOfInput() {
        final SAMFileHeader header = headerWithSample("s1");
        final int cap = 150;
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(cap, 300, header, new Random(1));
        // First input: a 1,000x block. Second input: overlapping, earlier-starting reads at 100x (under the cap),
        // as happens when one downsampler serves consecutive padded shards.
        final List<GATKRead> deep = pileOfReads(header, 500, 800, 10);
        final List<GATKRead> shallow = pileOfReads(header, 400, 900, 1);

        final List<GATKRead> keptDeep = runDownsampler(downsampler, deep);
        final List<GATKRead> keptShallow = runDownsampler(downsampler, shallow);

        Assert.assertTrue(keptDeep.size() < deep.size());
        // Nothing from the first input may count against the second: it is under the cap everywhere, so all of it survives.
        Assert.assertEquals(keptShallow, shallow);
        // Calling end-of-input again must be harmless.
        downsampler.signalEndOfInput();
        Assert.assertFalse(downsampler.hasFinalizedItems());
    }

    @Test
    public void testReadsReachingBeyondTheTrackedSpanAreAlwaysKept() {
        final SAMFileHeader header = headerWithSample("s1");
        final int windowSize = 100;
        final List<GATKRead> reads = new ArrayList<>(pileOfReads(header, 1, 100, 10)); // 1,000x, far over the cap
        // A spliced-style read whose span (100M then a 5 kb skip then 100M) reaches far past the tracked region.
        final GATKRead longSpan = ArtificialReadUtils.createArtificialRead(header, "spliced", 0, 50,
                new byte[200], new byte[200], "100M5000N100M");
        longSpan.setReadGroup(header.getReadGroups().get(0).getId());
        reads.add(longSpan);
        reads.sort(new ReadCoordinateComparator(header));

        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(20, windowSize, header, new Random(2)), reads);

        Assert.assertTrue(names(out).contains("spliced"), "a read reaching beyond the tracked span must be kept");
        Assert.assertTrue(out.size() < reads.size());
        assertCoordinateSorted(out, header);
    }

    @Test
    public void testSignalNoMoreReadsBeforeFinalizesTheCurrentWindow() {
        final SAMFileHeader header = headerWithSample("s1");
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(10, 100, header, new Random(1));
        downsampler.submit(read(header, "a", 1));
        downsampler.submit(read(header, "b", 50));
        Assert.assertTrue(downsampler.hasPendingItems());
        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertEquals(downsampler.size(), 2);
        Assert.assertEquals(downsampler.peekPending().getName(), "a");

        // A read still inside the window must not finalize it; one beyond the window must.
        downsampler.signalNoMoreReadsBefore(read(header, "inside", 90));
        Assert.assertTrue(downsampler.hasPendingItems());
        downsampler.signalNoMoreReadsBefore(read(header, "beyond", 101));
        Assert.assertFalse(downsampler.hasPendingItems());
        Assert.assertTrue(downsampler.hasFinalizedItems());
        Assert.assertEquals(downsampler.peekFinalized().getName(), "a");
        Assert.assertEquals(downsampler.consumeFinalizedItems().size(), 2);
        Assert.assertEquals(downsampler.size(), 0);
    }

    @Test
    public void testClearItemsDropsEverything() {
        final SAMFileHeader header = headerWithSample("s1");
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(10, 100, header, new Random(1));
        downsampler.submit(read(header, "a", 1));
        downsampler.submit(read(header, "b", 500));   // finalizes "a"
        Assert.assertTrue(downsampler.hasFinalizedItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.clearItems();

        Assert.assertFalse(downsampler.hasFinalizedItems());
        Assert.assertFalse(downsampler.hasPendingItems());
        Assert.assertEquals(downsampler.size(), 0);
        Assert.assertNull(downsampler.peekFinalized());
        Assert.assertNull(downsampler.peekPending());
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testUnsortedInputIsRejected() {
        final SAMFileHeader header = headerWithSample("s1");
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(10, 300, header, new Random(1));
        downsampler.submit(read(header, "a", 100));
        downsampler.submit(read(header, "b", 50));
    }

    @DataProvider(name = "invalidArguments")
    public Object[][] invalidArguments() {
        return new Object[][] {
                { 0, 300 },
                { -1, 300 },
                { 10, 0 },
                { 10, -5 },
        };
    }

    @Test(dataProvider = "invalidArguments", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidArgumentsAreRejected(final int cap, final int windowSize) {
        new MaxDepthDownsampler(cap, windowSize, ArtificialReadUtils.createArtificialSamHeader(), new Random(1));
    }
}
