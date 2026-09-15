package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

public final class MaxDepthDownsamplerUnitTest extends GATKBaseTest {

    private static final int READ_LENGTH = 100;

    private SAMFileHeader headerWithSample(final String sample) {
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg-" + sample);
        rg.setSample(sample);
        return ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(rg);
    }

    private GATKRead read(final SAMFileHeader header, final String name, final int start) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, start, READ_LENGTH);
        if (!header.getReadGroups().isEmpty()) {
            read.setReadGroup(header.getReadGroups().get(0).getId());
        }
        return read;
    }

    /** Coordinate-sorted reads at each start in [firstStart, lastStart], depthPerStart of them per start. */
    private List<GATKRead> pileOfReads(final SAMFileHeader header, final int firstStart, final int lastStart, final int depthPerStart) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int start = firstStart; start <= lastStart; start++) {
            for (int i = 0; i < depthPerStart; i++) {
                reads.add(read(header, "r" + start + "_" + i, start));
            }
        }
        return reads;
    }

    private List<GATKRead> runDownsampler(final MaxDepthDownsampler downsampler, final List<GATKRead> reads) {
        final List<GATKRead> out = new ArrayList<>();
        for (final GATKRead read : reads) {
            downsampler.submit(read);
            out.addAll(downsampler.consumeFinalizedItems());
        }
        downsampler.signalEndOfInput();
        out.addAll(downsampler.consumeFinalizedItems());
        return out;
    }

    private Map<Integer, Integer> depthByPosition(final List<GATKRead> reads) {
        final Map<Integer, Integer> depth = new TreeMap<>();
        for (final GATKRead read : reads) {
            for (int p = read.getStart(); p <= read.getEnd(); p++) {
                depth.merge(p, 1, Integer::sum);
            }
        }
        return depth;
    }

    private void assertSorted(final List<GATKRead> reads, final SAMFileHeader header) {
        for (int i = 1; i < reads.size(); i++) {
            Assert.assertTrue(ReadCoordinateComparator.compareCoordinates(reads.get(i - 1), reads.get(i), header) <= 0,
                    "output not coordinate sorted at index " + i);
        }
    }

    @Test
    public void reads_below_the_cap_pass_through_unchanged() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = pileOfReads(header, 1, 500, 1); // depth peaks at READ_LENGTH = 100
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(100, 300, header, new Random(1));
        final List<GATKRead> out = runDownsampler(downsampler, reads);
        Assert.assertEquals(out, reads);
        Assert.assertEquals(downsampler.getNumberOfDiscardedItems(), 0);
    }

    @Test
    public void depth_over_the_cap_is_reduced_but_never_below_the_cap() {
        final SAMFileHeader header = headerWithSample("s1");
        final int perStart = 20; // original depth 2,000x in the interior
        final List<GATKRead> reads = pileOfReads(header, 1, 1000, perStart);
        final int cap = 150;
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(cap, 300, header, new Random(7));
        final List<GATKRead> out = runDownsampler(downsampler, reads);

        final Map<Integer, Integer> original = depthByPosition(reads);
        final Map<Integer, Integer> kept = depthByPosition(out);
        for (final Map.Entry<Integer, Integer> e : original.entrySet()) {
            final int floor = Math.min(e.getValue(), cap);
            final int actual = kept.getOrDefault(e.getKey(), 0);
            Assert.assertTrue(actual >= floor, "position " + e.getKey() + " dropped to " + actual + " below " + floor);
        }
        // The cap should have bitten hard: kept depth must be far below the original in the interior.
        Assert.assertTrue(kept.get(500) < 3 * cap, "kept depth at 500 was " + kept.get(500));
        Assert.assertTrue(out.size() < reads.size() / 4, "expected most reads discarded, kept " + out.size() + " of " + reads.size());
        Assert.assertEquals(out.size() + downsampler.getNumberOfDiscardedItems(), reads.size());
        assertSorted(out, header);
    }

    @Test
    public void positions_under_the_cap_keep_every_read_even_next_to_deep_positions() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = new ArrayList<>();
        reads.addAll(pileOfReads(header, 1, 200, 1));      // shallow flank
        reads.addAll(pileOfReads(header, 201, 400, 10));   // deep block
        reads.addAll(pileOfReads(header, 401, 600, 1));    // shallow flank
        reads.sort(new ReadCoordinateComparator(header));
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(120, 300, header, new Random(3));
        final List<GATKRead> out = runDownsampler(downsampler, reads);
        final Set<String> keptNames = out.stream().map(GATKRead::getName).collect(Collectors.toSet());
        // Every read that touches a position whose original depth is under the cap must survive.
        final Map<Integer, Integer> original = depthByPosition(reads);
        for (final GATKRead read : reads) {
            boolean touchesShallow = false;
            for (int p = read.getStart(); p <= read.getEnd(); p++) {
                if (original.get(p) < 120) { touchesShallow = true; break; }
            }
            if (touchesShallow) {
                Assert.assertTrue(keptNames.contains(read.getName()), "read " + read.getName() + " touching a shallow position was discarded");
            }
        }
        Assert.assertTrue(downsampler.getNumberOfDiscardedItems() > 0);
    }

    @Test
    public void different_seeds_keep_different_subsets() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = pileOfReads(header, 1, 600, 10);
        final Set<String> a = runDownsampler(new MaxDepthDownsampler(100, 300, header, new Random(1)), reads).stream().map(GATKRead::getName).collect(Collectors.toSet());
        final Set<String> b = runDownsampler(new MaxDepthDownsampler(100, 300, header, new Random(2)), reads).stream().map(GATKRead::getName).collect(Collectors.toSet());
        Assert.assertNotEquals(a, b);
        final Set<String> a2 = runDownsampler(new MaxDepthDownsampler(100, 300, header, new Random(1)), reads).stream().map(GATKRead::getName).collect(Collectors.toSet());
        Assert.assertEquals(a, a2, "same seed must give the same subset");
    }

    @Test
    public void arrival_order_mode_is_deterministic_and_keeps_earliest_reads() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = pileOfReads(header, 1, 300, 10);
        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(50, 300, header, null), reads);
        final List<GATKRead> out2 = runDownsampler(new MaxDepthDownsampler(50, 300, header, null), reads);
        Assert.assertEquals(out, out2);
        // In arrival order the first 50 reads (5 starts x 10) are always kept.
        for (int i = 0; i < 50; i++) {
            Assert.assertEquals(out.get(i).getName(), reads.get(i).getName());
        }
    }

    @Test
    public void depth_is_capped_per_sample() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1"); rg1.setSample("s1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2"); rg2.setSample("s2");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.addReadGroup(rg1);
        header.addReadGroup(rg2);
        final List<GATKRead> reads = new ArrayList<>();
        for (int start = 1; start <= 300; start++) {
            for (int i = 0; i < 4; i++) {
                final GATKRead r = ArtificialReadUtils.createArtificialRead(header, "s1_" + start + "_" + i, 0, start, READ_LENGTH);
                r.setReadGroup("rg1"); reads.add(r);
            }
            final GATKRead r = ArtificialReadUtils.createArtificialRead(header, "s2_" + start, 0, start, READ_LENGTH);
            r.setReadGroup("rg2"); reads.add(r);
        }
        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(150, 300, header, new Random(5)), reads);
        final List<GATKRead> s2 = out.stream().filter(r -> r.getReadGroup().equals("rg2")).collect(Collectors.toList());
        // Sample 2 peaks at 100x, under the cap, so all of its reads survive regardless of sample 1's depth.
        Assert.assertEquals(s2.size(), 300);
        Assert.assertTrue(out.size() < reads.size());
    }

    @Test
    public void reads_without_an_assigned_position_pass_through() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = new ArrayList<>(pileOfReads(header, 1, 50, 1));
        final GATKRead unmapped = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30});
        reads.add(unmapped);
        final List<GATKRead> out = runDownsampler(new MaxDepthDownsampler(10, 300, header, new Random(1)), reads);
        Assert.assertTrue(out.contains(unmapped));
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void unsorted_input_is_rejected() {
        final SAMFileHeader header = headerWithSample("s1");
        final MaxDepthDownsampler downsampler = new MaxDepthDownsampler(10, 300, header, new Random(1));
        downsampler.submit(read(header, "a", 100));
        downsampler.submit(read(header, "b", 50));
    }

    @Test
    public void chained_downsampler_applies_both_stages() {
        final SAMFileHeader header = headerWithSample("s1");
        final List<GATKRead> reads = pileOfReads(header, 1, 300, 10);
        final ChainedReadsDownsampler chained = new ChainedReadsDownsampler(
                new PositionalDownsampler(5, header, true),
                new MaxDepthDownsampler(50, 300, header, new Random(1)));
        final List<GATKRead> out = new ArrayList<>();
        for (final GATKRead r : reads) {
            chained.submit(r);
            out.addAll(chained.consumeFinalizedItems());
        }
        chained.signalEndOfInput();
        out.addAll(chained.consumeFinalizedItems());
        final Map<Integer, Integer> kept = depthByPosition(out);
        Assert.assertTrue(kept.get(150) >= 50 && kept.get(150) < 200, "depth at 150 was " + kept.get(150));
        Assert.assertEquals(chained.getNumberOfDiscardedItems() + out.size(), reads.size());
        assertSorted(out, header);
        Assert.assertEquals(new HashSet<>(out).size(), out.size());
    }
}
