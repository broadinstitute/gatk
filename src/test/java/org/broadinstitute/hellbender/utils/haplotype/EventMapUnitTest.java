package org.broadinstitute.hellbender.utils.haplotype;


import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.UnvalidatingGenomeLoc;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class EventMapUnitTest extends GATKBaseTest {
    private final static String CHR = "20";
    private final static String NAME = "foo";

    @DataProvider(name = "MNPTest")
    public Object[][] makeMNPTest() {
        List<Object[]> tests = new ArrayList<>();
        final List<Integer> justZero = Collections.singletonList(0);

        tests.add(new Object[]{"TTTGGGAAA", "TTTCCCAAA", "3M3X3M", Arrays.asList(1, 2, 3, 5, 10),
                Arrays.asList(Arrays.asList("GGG", "CCC"))});
        tests.add(new Object[]{"TTTGGGAAA", "TTTCCCAAA", "3M3X3M", justZero ,
                Arrays.asList(Arrays.asList("G", "C"), Arrays.asList("G", "C"), Arrays.asList("G", "C"))});
        tests.add(new Object[]{"TTTGGGAAA", "TTTCCCAAA", "9M", Arrays.asList(1, 2, 3, 5, 10),
                Arrays.asList(Arrays.asList("GGG", "CCC"))});
        tests.add(new Object[]{"TTTGGGAAA", "TTTCCCAAA", "9M", justZero,
                Arrays.asList(Arrays.asList("G", "C"), Arrays.asList("G", "C"), Arrays.asList("G", "C"))});
        tests.add(new Object[]{"TTTTTTTTT", "ATATATATA", "9M", Arrays.asList(2),
                Arrays.asList(Arrays.asList("TTTTTTTTT", "ATATATATA"))});
        tests.add(new Object[]{"ACGT", "CGTA", "4M", Arrays.asList(1, 2, 3, 5, 10),
                Arrays.asList(Arrays.asList("ACGT", "CGTA"))});
        tests.add(new Object[]{"ACGT", "CGTA", "4M", justZero,
                Arrays.asList(Arrays.asList("A", "C"), Arrays.asList("C", "G"), Arrays.asList("G", "T"), Arrays.asList("T", "A"))});
        tests.add(new Object[]{"ACTTGC", "CATTCG", "6M", Arrays.asList(1, 2),
                Arrays.asList(Arrays.asList("AC", "CA"), Arrays.asList("GC", "CG"))});
        tests.add(new Object[]{"ACTTGC", "CATTCG", "6M", Arrays.asList(3,5,10),
                Arrays.asList(Arrays.asList("ACTTGC", "CATTCG"))});
        tests.add(new Object[]{"ACTTGC", "CATTCG", "6M", justZero,
                Arrays.asList(Arrays.asList("A", "C"), Arrays.asList("C", "A"), Arrays.asList("G", "C"), Arrays.asList("C", "G"))});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MNPTest")
    public void testMNPs(final String refBases, final String haplotypeBases, final String cigar, final List<Integer> maxMnpDistance, final List<List<String>> expectedAlleles) {
        final Haplotype hap = new Haplotype(haplotypeBases.getBytes(), false, 0, TextCigarCodec.decode(cigar));
        final GenomeLoc loc = new UnvalidatingGenomeLoc(CHR, 0, 1, refBases.length());
        for (final int maxDist : maxMnpDistance) {
            final EventMap events = new EventMap(hap, refBases.getBytes(), loc, NAME, maxDist);
            Assert.assertEquals(events.getNumberOfEvents(), expectedAlleles.size());
            final List<VariantContext> foundAlleles = new ArrayList<>(events.getVariantContexts());
            for (int i = 0; i < events.getNumberOfEvents(); i++) {
                final VariantContext actual = foundAlleles.get(i);
                Assert.assertEquals(actual.getReference().getDisplayString(), expectedAlleles.get(i).get(0));
                Assert.assertEquals(actual.getAlternateAllele(0).getDisplayString(), expectedAlleles.get(i).get(1));
            }
        }
    }

    @DataProvider(name = "MakeBlockData")
    public Object[][] makeMakeBlockData() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{Arrays.asList("A", "G"), Arrays.asList("AGT", "A"), Arrays.asList("AGT", "G")});
        tests.add(new Object[]{Arrays.asList("A", "G"), Arrays.asList("A", "AGT"), Arrays.asList("A", "GGT")});

        tests.add(new Object[]{Arrays.asList("AC", "A"), Arrays.asList("A", "AGT"), Arrays.asList("AC", "AGT")});
        tests.add(new Object[]{Arrays.asList("ACGTA", "A"), Arrays.asList("A", "AG"), Arrays.asList("ACGTA", "AG")});
        tests.add(new Object[]{Arrays.asList("AC", "A"), Arrays.asList("A", "AGCGT"), Arrays.asList("AC", "AGCGT")});
        tests.add(new Object[]{Arrays.asList("A", "ACGTA"), Arrays.asList("AG", "A"), Arrays.asList("AG", "ACGTA")});
        tests.add(new Object[]{Arrays.asList("A", "AC"), Arrays.asList("AGCGT", "A"), Arrays.asList("AGCGT", "AC")});

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "MakeBlockData")
    public void testGetNeighborhood(final List<String> firstAlleles, final List<String> secondAlleles, final List<String> expectedAlleles) {
        final VariantContext vc1 = GATKVariantContextUtils.makeFromAlleles("x", "20", 10, firstAlleles);
        final VariantContext vc2 = GATKVariantContextUtils.makeFromAlleles("x", "20", 10, secondAlleles);
        final VariantContext expected = GATKVariantContextUtils.makeFromAlleles("x", "20", 10, expectedAlleles);

        final EventMap eventMap = new EventMap(Collections.<VariantContext>emptyList());
        final VariantContext block = eventMap.makeBlock(vc1, vc2);

        Assert.assertEquals(block.getStart(), expected.getStart());
        Assert.assertEquals(block.getAlleles(), expected.getAlleles());
    }
}
