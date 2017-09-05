package org.broadinstitute.hellbender.utils.collections;

import com.google.common.collect.Lists;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class IntervalsSkipListUnitTest extends GATKBaseTest {

    @DataProvider(name="intervals")
    public Object[][] intervals(){
        ArrayList<Locatable> input = Lists.newArrayList(
                new SimpleInterval("1",  10, 100),
                new SimpleInterval("2", 200, 300)
        );
        ArrayList<Locatable> empty = new ArrayList<>();
        ArrayList<Locatable> contig1 = Lists.newArrayList(
                new SimpleInterval("1",  10, 100)
        );
        ArrayList<Locatable> contig2 = Lists.newArrayList(
                new SimpleInterval("2", 200, 300)
                );

        // returns input, query range, expected SimpleIntervals
        return new Object[][]{
                // we already test elsewhere that it works within a contig, so here we just have to make sure that
                // it picks the correct contig and can deal with not-yet-mentioned contigs.
                new Object[]{input, new SimpleInterval("1", 100, 200), contig1},
                new Object[]{input, new SimpleInterval("1", 1, 5), empty},
                new Object[]{input, new SimpleInterval("2", 100, 200), contig2},
                new Object[]{input, new SimpleInterval("3", 100, 200), empty},
        };
    }

    @DataProvider(name = "intervalSets")
    public Object[][] intervalSets() {
        final List<Set<? extends Locatable>> result = new ArrayList<>();
        final Random rdn = new Random(13);
        result.add(Collections.emptySet());
        result.add(Collections.singleton(new SimpleInterval("seq0", 100, 200)));
        for (int i = 0; i < 1_000; i++) {
            final Set<SimpleInterval> intervals = new HashSet<>();
            final int contigCount = rdn.nextInt(3) + 1;
            final int[] contigSizes = IntStream.range(0, contigCount).map(c -> rdn.nextInt(1000) + 10).toArray();
            final int[] accumulateSizes = new int[contigCount];
            for (int j = 1; j < contigCount; j++) {
                accumulateSizes[j] = accumulateSizes[j - 1] + contigSizes[j - 1];
            }
            final double depth = rdn.nextDouble() * 2;
            long totalBases = 0;
            while (totalBases < accumulateSizes[contigCount - 1] * depth) {
                int start = rdn.nextInt(accumulateSizes[contigCount - 1]);
                int contigIndex = 0;
                for (int j = 0; j < contigCount - 1; j++, contigIndex++) {
                    if (start < accumulateSizes[j + 1]) {
                        break;
                    }
                    start -= contigSizes[j];
                }
                final int length = rdn.nextInt(50)  + 1;
                final int end = Math.min(start + length, contigSizes[contigIndex]);
                final SimpleInterval si = new SimpleInterval("seq" + contigIndex, start + 1, end);
                if (intervals.add(si)) {
                    totalBases += length;
                }
            }
            result.add(intervals);
        }
        return result.stream().map(l -> new Object[] { l }).toArray(Object[][]::new);
    }

    @Test(dataProvider = "intervals")
    public void testOverlap(ArrayList<Locatable> input, SimpleInterval query, ArrayList<Locatable> expected) throws Exception {
        IntervalsSkipList<Locatable> ints = new IntervalsSkipList<>(input);
        List<Locatable> actual = ints.getOverlapping(query);
        Assert.assertEquals(
                actual,
                expected
        );
    }

    @Test(dataProvider = "intervalSets")
    public void testSize(final Set<? extends Locatable> input) throws Exception {
        final IntervalsSkipList<? extends Locatable> ints = new IntervalsSkipList<>(input);
        Assert.assertEquals(ints.size(), input.stream().distinct().count());
    }

    @Test(dataProvider= "intervalSets")
    public void testContigs(final Set<? extends Locatable> input) throws Exception {
        final IntervalsSkipList<? extends Locatable> ints = new IntervalsSkipList<>(input);
        final List<String> expectedContigs = input.stream().map(Locatable::getContig)
                .collect(Collectors.toList());
        final List<String> actualContigs = ints.contigs();
        Assert.assertEquals(actualContigs, expectedContigs);
    }

    @Test(dataProvider = "intervalSets", dependsOnMethods = "testContigs")
    public void testContigIterators(final Set<? extends Locatable> input) throws Exception {
        final IntervalsSkipList<? extends Locatable> ints = new IntervalsSkipList<>(input);
        for (final String contig : ints.contigs()) {
            final Iterator<? extends Locatable> contigLocatables = ints.contigIterator(contig);
            Assert.assertNotNull(contigLocatables);
            final List<Locatable> contigActual = Utils.stream(contigLocatables)
                    .collect(Collectors.toList());
            final List<Locatable> contigExpected = input.stream()
                    .filter(l -> contig.equals(l.getContig()))
                    .sorted(Comparator.comparingInt(Locatable::getStart)
                                      .thenComparingInt(Locatable::getEnd))
                    .collect(Collectors.toList());
            Assert.assertEquals(contigActual, contigExpected);
        }
        // check on an non-existent contig:
        Assert.assertFalse(ints.contigIterator("non-existent-contig").hasNext());
    }


}