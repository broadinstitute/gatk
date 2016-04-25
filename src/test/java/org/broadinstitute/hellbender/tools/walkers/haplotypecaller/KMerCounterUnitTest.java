package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public final class KMerCounterUnitTest extends BaseTest {
    @Test
	public void testMyData() {
        final KMerCounter counter = new KMerCounter(3);

        Assert.assertNotNull(counter.toString());

        counter.addKmers(
			 "ATG", "ATG", "ATG", "ATG",
			 "ACC", "ACC", "ACC",
			 "AAA", "AAA",
			 "CTG",
			 "NNA",
                "CCC"
			 );

        testCounting(counter, "ATG", 4);
        testCounting(counter, "ACC", 3);
        testCounting(counter, "AAA", 2);
        testCounting(counter, "CTG", 1);
        testCounting(counter, "NNA", 1);
        testCounting(counter, "CCC", 1);
        testCounting(counter, "NNN", 0);
        testCounting(counter, "NNC", 0);

        Assert.assertNotNull(counter.toString());

        assertCounts(counter, 5);
        assertCounts(counter, 4, "ATG");
        assertCounts(counter, 3, "ATG", "ACC");
        assertCounts(counter, 2, "ATG", "ACC", "AAA");
        assertCounts(counter, 1, "ATG", "ACC", "AAA", "CTG", "NNA", "CCC");

        counter.clear();
        assertCounts(counter, 0);
    }

    private void assertCounts(final KMerCounter counter, final int minCount, final String... expecteds) {
        final Set<Kmer> expected = new HashSet<>();
        for ( final String one : expecteds ) expected.add(new Kmer(one));
        Assert.assertEquals(new HashSet<>(counter.getKmersWithCountsAtLeast(minCount)), expected);
    }

    private void testCounting(final KMerCounter counter, final String in, final int expectedCount) {
        Assert.assertEquals(counter.getKmerCount(new Kmer(in)), expectedCount);
    }

    @Test
    public void testCountedKMers() throws Exception {
        final KMerCounter counter = new KMerCounter(3);

        final String kmer = "ATG";
        counter.addKmers(
                kmer
        );
        final Collection<KMerCounter.CountedKmer> countedKmers = counter.getCountedKmers();
        Assert.assertEquals(countedKmers.size(), 1);
        KMerCounter.CountedKmer ckmer = countedKmers.iterator().next();
        Assert.assertNotNull(ckmer.toString());
        Assert.assertEquals(ckmer.getKmer().bases(), kmer.getBytes());
        Assert.assertEquals(ckmer.getCount(), 1);
    }

    @Test
    public void testCountedKMersCompare() throws Exception {
        final KMerCounter counter = new KMerCounter(3);

        final String kmer1 = "ATG";
        final String kmer2 = "CTG";
        counter.addKmers(
                kmer1,
                kmer2, kmer2
        );
        final Collection<KMerCounter.CountedKmer> countedKmers = counter.getCountedKmers();
        Assert.assertEquals(countedKmers.size(), 2);
        final List<KMerCounter.CountedKmer> list = new ArrayList<>(countedKmers);
        list.sort((ck1, ck2) -> ck1.compareTo(ck2));
        Assert.assertEquals(list.get(0).getKmer().bases(), kmer2.getBytes());
        Assert.assertEquals(list.get(1).getKmer().bases(), kmer1.getBytes());
    }
}
