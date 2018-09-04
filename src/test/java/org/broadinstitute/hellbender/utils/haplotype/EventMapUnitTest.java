package org.broadinstitute.hellbender.utils.haplotype;


import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.GenomeLoc;
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

    @DataProvider(name = "MyDataProvider")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final List<String> SNP_ALLELES = Arrays.asList("A", "C");
        final List<String> INS_ALLELES = Arrays.asList("A", "ACGTGA");
        final List<String> DEL_ALLELES = Arrays.asList("ACGTA", "C");
        final List<List<String>> allAlleles = Arrays.asList(SNP_ALLELES, INS_ALLELES, DEL_ALLELES);
        for ( final int leftNotClump : Arrays.asList(-1, 3) ) {
            for ( final int middleNotClump : Arrays.asList(-1, 10, 500) ) {
                for ( final int rightNotClump : Arrays.asList(-1, 1000) ) {
                    for ( final int nClumped : Arrays.asList(3, 4) ) {
                        for ( final List<List<String>> alleles : Utils.makePermutations(allAlleles, nClumped, true)) {
                            final List<VariantContext> allVCS = new LinkedList<>();

                            if ( leftNotClump != -1 ) allVCS.add(GATKVariantContextUtils.makeFromAlleles(NAME, CHR, leftNotClump, SNP_ALLELES));
                            if ( middleNotClump != -1 ) allVCS.add(GATKVariantContextUtils.makeFromAlleles(NAME, CHR, middleNotClump, SNP_ALLELES));
                            if ( rightNotClump != -1 ) allVCS.add(GATKVariantContextUtils.makeFromAlleles(NAME, CHR, rightNotClump, SNP_ALLELES));

                            int clumpStart = 50;
                            final List<VariantContext> vcs = new LinkedList<>();
                            for ( final List<String> myAlleles : alleles ) {
                                final VariantContext vc = GATKVariantContextUtils.makeFromAlleles(NAME, CHR, clumpStart, myAlleles);
                                clumpStart = vc.getEnd() + 3;
                                vcs.add(vc);
                            }

                            tests.add(new Object[]{new EventMap(new LinkedList<>(allVCS)), Collections.emptyList()});
                            allVCS.addAll(vcs);
                            tests.add(new Object[]{new EventMap(allVCS), vcs});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "MyDataProvider")
    public void testGetNeighborhood(final EventMap eventMap, final List<VariantContext> expectedNeighbors) {
        final VariantContext leftOfNeighors = expectedNeighbors.isEmpty() ? null : expectedNeighbors.get(0);

        for ( final VariantContext vc : eventMap.getVariantContexts() ) {
            final List<VariantContext> n = eventMap.getNeighborhood(vc, 5);
            if ( leftOfNeighors == vc )
                Assert.assertEquals(n, expectedNeighbors);
            else if ( ! expectedNeighbors.contains(vc) )
                Assert.assertEquals(n, Collections.singletonList(vc), "Should only contain the original vc but " + n);
        }
    }

    @DataProvider(name = "BlockSubstitutionsData")
    public Object[][] makeBlockSubstitutionsData() {
        List<Object[]> tests = new ArrayList<>();

        for ( int size = EventMap.MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION; size < 10; size++ ) {
            final String ref = StringUtils.repeat("A", size);
            final String alt = StringUtils.repeat("C", size);
            tests.add(new Object[]{ref, alt, size + "M", GATKVariantContextUtils.makeFromAlleles(NAME, CHR, 1, Arrays.asList(ref, alt))});
        }

        tests.add(new Object[]{"AAAAAA", "GAGAGA", "6M", GATKVariantContextUtils.makeFromAlleles(NAME, CHR, 1, Arrays.asList("AAAAA", "GAGAG"))});
        tests.add(new Object[]{"AAAAAA", "GAGAGG", "6M", GATKVariantContextUtils.makeFromAlleles(NAME, CHR, 1, Arrays.asList("AAAAAA", "GAGAGG"))});

        for ( int len = 0; len < 10; len++ ) {
            final String s = len == 0 ? "" : StringUtils.repeat("A", len);
            tests.add(new Object[]{s + "AACCCCAA", s + "GAAG", len + 2 + "M4D2M", GATKVariantContextUtils.makeFromAlleles(NAME, CHR, 1 + len,   Arrays.asList("AACCCCAA", "GAAG"))});
            tests.add(new Object[]{s + "AAAA", s + "GACCCCAG", len + 2 + "M4I2M", GATKVariantContextUtils.makeFromAlleles(NAME, CHR, 1 + len, Arrays.asList("AAAA", "GACCCCAG"))});

            tests.add(new Object[]{"AACCCCAA" + s, "GAAG" + s, "2M4D" + (len + 2) + "M", GATKVariantContextUtils.makeFromAlleles(NAME, CHR, 1,   Arrays.asList("AACCCCAA", "GAAG"))});
            tests.add(new Object[]{"AAAA" + s, "GACCCCAG" + s, "2M4I" + (len + 2) + "M", GATKVariantContextUtils.makeFromAlleles(NAME, CHR, 1, Arrays.asList("AAAA", "GACCCCAG"))});
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "BlockSubstitutionsData")
    public void testBlockSubstitutionsData(final String refBases, final String haplotypeBases, final String cigar, final VariantContext expectedBlock) {
        final Haplotype hap = new Haplotype(haplotypeBases.getBytes(), false, 0, TextCigarCodec.decode(cigar));
        final GenomeLoc loc = new UnvalidatingGenomeLoc(CHR, 0, 1, refBases.length());
        final EventMap ee = new EventMap(hap, refBases.getBytes(), loc, NAME, 1);
        ee.replaceClumpedEventsWithBlockSubstitutions();
        Assert.assertEquals(ee.getNumberOfEvents(), 1);
        final VariantContext actual = ee.getVariantContexts().iterator().next();
        Assert.assertTrue(GATKVariantContextUtils.equalSites(actual, expectedBlock), "Failed with " + actual);
    }

    @DataProvider(name = "AdjacentSNPIndelTest")
    public Object[][] makeAdjacentSNPIndelTest() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{"TT", "GCT", "1M1I1M", Arrays.asList(Arrays.asList("T", "GC"))});
        tests.add(new Object[]{"GCT", "TT", "1M1D1M", Arrays.asList(Arrays.asList("GC", "T"))});
        tests.add(new Object[]{"TT", "GCCT", "1M2I1M", Arrays.asList(Arrays.asList("T", "GCC"))});
        tests.add(new Object[]{"GCCT", "TT", "1M2D1M", Arrays.asList(Arrays.asList("GCC", "T"))});
        tests.add(new Object[]{"AAGCCT", "AATT", "3M2D1M", Arrays.asList(Arrays.asList("GCC", "T"))});
        tests.add(new Object[]{"AAGCCT", "GATT", "3M2D1M", Arrays.asList(Arrays.asList("A", "G"), Arrays.asList("GCC", "T"))});
        tests.add(new Object[]{"AAAAA", "AGACA", "5M", Arrays.asList(Arrays.asList("A", "G"), Arrays.asList("A", "C"))});

        return tests.toArray(new Object[][]{});
    }

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

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "AdjacentSNPIndelTest")
    public void testAdjacentSNPIndelTest(final String refBases, final String haplotypeBases, final String cigar, final List<List<String>> expectedAlleles) {
        final Haplotype hap = new Haplotype(haplotypeBases.getBytes(), false, 0, TextCigarCodec.decode(cigar));
        final GenomeLoc loc = new UnvalidatingGenomeLoc(CHR, 0, 1, refBases.length());
        final EventMap ee = new EventMap(hap, refBases.getBytes(), loc, NAME, 1);
        ee.replaceClumpedEventsWithBlockSubstitutions();
        Assert.assertEquals(ee.getNumberOfEvents(), expectedAlleles.size());
        final List<VariantContext> actuals = new ArrayList<>(ee.getVariantContexts());
        for ( int i = 0; i < ee.getNumberOfEvents(); i++ ) {
            final VariantContext actual = actuals.get(i);
            Assert.assertEquals(actual.getReference().getDisplayString(), expectedAlleles.get(i).get(0));
            Assert.assertEquals(actual.getAlternateAllele(0).getDisplayString(), expectedAlleles.get(i).get(1));
        }
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
