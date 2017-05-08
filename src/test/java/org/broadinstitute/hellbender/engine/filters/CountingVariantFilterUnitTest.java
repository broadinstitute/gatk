package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class CountingVariantFilterUnitTest {

    static final VariantContext goodVariant = new VariantContextBuilder("Zuul", "1", 2, 2, Collections.singletonList(Allele.create("A", true))).make();
    static final VariantContext endBad = new VariantContextBuilder("Peter", "1", 2, 20, Collections.singletonList(Allele.create("TTTTTTTTTTTTTTTTTTT", true))).make();
    static final VariantContext startBad = new VariantContextBuilder("Ray", "1", 1, 2, Collections.singletonList(Allele.create("AA", true))).make();
    static final VariantContext bothBad = new VariantContextBuilder("Egon", "1", 1, 20, Collections.singletonList(Allele.create("TTTTTTTTTTTTTTTTTTTT", true))).make();

    static final VariantFilter startOk = new VariantFilter() {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final VariantContext variant){return variant.getStart() > 1;}
    };
    static final VariantFilter endOk = new VariantFilter() {
        private static final long serialVersionUID = 1L;
        @Override public boolean test(final VariantContext variant){return variant.getEnd() <= 10;}
    };

    // Helper to verify post-filtering filter state
    private void verifyFilterState(CountingVariantFilter vf, boolean expected) {
        long count = vf.getFilteredCount();
        String vfSummary = vf.getSummaryLine();

        if (expected) {
            Assert.assertTrue(0 == count);
            Assert.assertEquals(0, vfSummary.indexOf("No variants filtered"));
        } else {
            Assert.assertTrue(1 == count);
            Assert.assertEquals(0, vfSummary.indexOf("1 variant(s) filtered"));
        }
    }

    @DataProvider(name = "variantsStartEnd")
    public Object[][] variantsStartEnd() {
        return new Object[][]{
                {goodVariant, true, true},
                {startBad, false, true},
                {endBad, true, false},
                {bothBad, false, false}
        };
    }

    @Test(dataProvider = "variantsStartEnd")
    public void testTest(VariantContext variant, boolean start, boolean end) {

        CountingVariantFilter startOkCounting = new CountingVariantFilter(startOk);
        Assert.assertEquals(startOkCounting.test(variant), start);
        verifyFilterState(startOkCounting, start);

        CountingVariantFilter endOkCounting = new CountingVariantFilter(endOk);
        Assert.assertEquals(endOkCounting.test(variant), end);
        verifyFilterState(endOkCounting, end);
    }

    @Test(dataProvider = "variantsStartEnd")
    public void testNegate(VariantContext variant, boolean start, boolean end) {
        CountingVariantFilter notStartOkCounting = new CountingVariantFilter(startOk).negate();
        Assert.assertEquals(notStartOkCounting.test(variant), !start);
        verifyFilterState(notStartOkCounting, !start);

        CountingVariantFilter notEndOkCounting = new CountingVariantFilter(endOk).negate();
        Assert.assertEquals(notEndOkCounting.test(variant), !end);
        verifyFilterState(notEndOkCounting, !end);
    }

    @DataProvider(name = "variantsAnd")
    public Object[][] variantsAnd() {
        return new Object[][]{
                {goodVariant, true},
                {startBad, false},
                {endBad, false},
                {bothBad, false}
        };
    }

    @Test(dataProvider = "variantsAnd")
    public void testAnd(VariantContext variant, boolean expected) {

        CountingVariantFilter startAndEndOk = new CountingVariantFilter(startOk).and(new CountingVariantFilter(endOk));
        Assert.assertEquals(startAndEndOk.test(variant), expected);
        verifyFilterState(startAndEndOk, expected);

        CountingVariantFilter endAndStartOk = new CountingVariantFilter(endOk).and(new CountingVariantFilter(startOk));
        Assert.assertEquals(endAndStartOk.test(variant), expected);
        verifyFilterState(endAndStartOk, expected);
    }

    @DataProvider(name = "variantsOr")
    public Object[][] variantsOr() {
        return new Object[][]{
                {goodVariant, true},
                {startBad, true},
                {endBad, true},
                {bothBad, false}
        };
    }

    @Test(dataProvider = "variantsOr")
    public void testOr(VariantContext variant, boolean expected) {

        CountingVariantFilter startOrEndOk = new CountingVariantFilter(startOk).or(new CountingVariantFilter(endOk));
        Assert.assertEquals(startOrEndOk.test(variant), expected);
        verifyFilterState(startOrEndOk, expected);

        CountingVariantFilter endOrStartOk = new CountingVariantFilter(endOk).or(new CountingVariantFilter(startOk));
        Assert.assertEquals(endOrStartOk.test(variant), expected);
        verifyFilterState(endOrStartOk, expected);
    }

    @DataProvider(name = "deeper")
    public Object[][] deeper() {
        return new Object[][]{
                {goodVariant, false},
                {startBad, true},
                {endBad, true},
                {bothBad, false}
        };
    }

    private CountingVariantFilter variantChecksOut() {
        return new CountingVariantFilter(startOk)
                .or(new CountingVariantFilter(endOk))
                .and(new CountingVariantFilter(new VariantFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final VariantContext variant){return !variant.getSource().equals("Zuul");}
                }));
    }

    @Test(dataProvider = "deeper")
    public void testDeeperChaining(VariantContext variant, boolean expected) {

        CountingVariantFilter variantCheckOutCounting = variantChecksOut();
        Assert.assertEquals(variantCheckOutCounting.test(variant), expected);
        verifyFilterState(variantCheckOutCounting, expected);

        variantCheckOutCounting = variantChecksOut().and(variantChecksOut());
        Assert.assertEquals(variantCheckOutCounting.test(variant), expected);
        verifyFilterState(variantCheckOutCounting, expected);

        variantCheckOutCounting = variantChecksOut().and(new CountingVariantFilter(
                new VariantFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final VariantContext variant){return false;}
                }));
        Assert.assertEquals(variantCheckOutCounting.test(variant), false);
        verifyFilterState(variantCheckOutCounting, false);

        variantCheckOutCounting = variantChecksOut().or(new CountingVariantFilter(
                new VariantFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final VariantContext variant){return true;}
                }));
        Assert.assertEquals(variantCheckOutCounting.test(variant), true);
        verifyFilterState(variantCheckOutCounting, true);
    }

    @DataProvider(name = "multipleRejection")
    public Object[][] multipleRejection() {
        return new Object[][] {
            {
                new VariantContext[] {goodVariant, goodVariant, goodVariant}, 0
            },
            {
                new VariantContext[] {goodVariant, goodVariant, bothBad }, 1
            },
            {
                new VariantContext[] { startBad, bothBad, goodVariant, startBad }, 3
            },
        };
    }

    @Test(dataProvider = "multipleRejection")
    public void testRootFilterCounts(VariantContext[] variants, int expectedRejectionCount) {

        CountingVariantFilter startOkCounting = new CountingVariantFilter(startOk);
        Arrays.asList(variants).stream().filter(startOkCounting).count();  // force the stream to be consumed
        Assert.assertEquals(startOkCounting.getFilteredCount(), expectedRejectionCount);
        startOkCounting.resetFilteredCount();
        Assert.assertEquals(startOkCounting.getFilteredCount(), 0);
    }

    @DataProvider(name = "subFilterCounts")
    public Object[][] subFilterCounts() {
        return new Object[][] {
                {
                        new VariantContext[]{goodVariant, startBad, bothBad, bothBad}, 1L, 2L, 1L
                },
                {
                        new VariantContext[]{goodVariant, goodVariant, goodVariant, bothBad }, 3L, 3L, 3L
                },
                {
                        new VariantContext[]{goodVariant, startBad, endBad, bothBad}, 2L, 3L, 2L
                },
        };
    }

    @Test(dataProvider = "subFilterCounts")
    public void testSubFilterCounts(VariantContext[] variants, long totalRejections, long startEndRejections, long nameRejections) {

        CountingVariantFilter badStart = new CountingVariantFilter(
                new VariantFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final VariantContext variant){return variant.getStart() <= 1;}
                }
        );
        CountingVariantFilter badEnd = new CountingVariantFilter(
                new VariantFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final VariantContext variant){return variant.getEnd() > 10;}
                }
        );
        CountingVariantFilter badStartAndEnd = badStart.and(badEnd);

        CountingVariantFilter isRay= new CountingVariantFilter(
                new VariantFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final VariantContext variant){return variant.getSource().equals("Ray");}
                }
        );
        CountingVariantFilter isEgon = new CountingVariantFilter(
                new VariantFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override public boolean test(final VariantContext variant){return variant.getSource().equals("Egon");}
                }
        );
        CountingVariantFilter isRayOrEgon = isRay.or(isEgon);

        CountingVariantFilter compoundFilter = badStartAndEnd.or(isRayOrEgon);

        Arrays.asList(variants).stream().filter(compoundFilter).count(); // force the stream to be consumed

        Assert.assertEquals(compoundFilter.getFilteredCount(), totalRejections);
        Assert.assertEquals(badStartAndEnd.getFilteredCount(), startEndRejections);
        Assert.assertEquals(isRayOrEgon.getFilteredCount(), nameRejections);

        // test if reset filtered count is correctly propagated
        compoundFilter.resetFilteredCount();
        Assert.assertEquals(compoundFilter.getFilteredCount(), 0);
        Assert.assertEquals(badStartAndEnd.getFilteredCount(), 0);
        Assert.assertEquals(isRayOrEgon.getFilteredCount(), 0);
        Assert.assertEquals(badStart.getFilteredCount(), 0);
        Assert.assertEquals(badEnd.getFilteredCount(), 0);
        Assert.assertEquals(isRay.getFilteredCount(), 0);
        Assert.assertEquals(isEgon.getFilteredCount(), 0);
    }

    @Test
    public void testFromListNull() {
        CountingVariantFilter vf = CountingVariantFilter.fromList(null);
        Assert.assertTrue(vf.delegateFilter == VariantFilterLibrary.ALLOW_ALL_VARIANTS);
    }
    @Test
    public void testFromListEmpty() {
        CountingVariantFilter vf = CountingVariantFilter.fromList(Collections.emptyList());
        Assert.assertTrue(vf.delegateFilter == VariantFilterLibrary.ALLOW_ALL_VARIANTS);
    }

    @Test
    public void testFromListSingle() {
        final VariantFilter snpFilter = VariantContext::isSNP;
        List<VariantFilter> filters = new ArrayList<>();
        filters.add(snpFilter);
        CountingVariantFilter vf = CountingVariantFilter.fromList(filters);
        Assert.assertTrue(vf.delegateFilter == snpFilter);
    }

    @Test
    public void testFromListMultiOrdered() {
        final VariantFilter snpFilter = VariantContext::isSNP;
        final VariantFilter biallelicFilter = VariantContext::isBiallelic;
        final VariantFilter symbolicFilter = VariantContext::isSymbolic;
        List<VariantFilter> filters = new ArrayList<>();
        filters.add(snpFilter);
        filters.add(biallelicFilter);
        filters.add(symbolicFilter);

        // Since we want to ensure that order of the input is honored, we need to test the
        // structure of the filter rather than the result
        CountingVariantFilter vf = CountingVariantFilter.fromList(filters);

        Assert.assertTrue(vf.getClass() == CountingVariantFilter.CountingAndVariantFilter.class);
        CountingVariantFilter.CountingAndVariantFilter andFilter = (CountingVariantFilter.CountingAndVariantFilter) vf;

        // lhs is a Counting and filter; rhs is a counting filter that delegates to symbolicFilter
        Assert.assertTrue(andFilter.lhs.getClass() == CountingVariantFilter.CountingAndVariantFilter.class);
        Assert.assertTrue(andFilter.rhs.delegateFilter == symbolicFilter);
        andFilter = (CountingVariantFilter.CountingAndVariantFilter) andFilter.lhs;

        // lhs is a Counting filter that delegates to snpFilter; rhs is a
        // counting filter that delegates to biallelicFilter
        Assert.assertTrue(andFilter.lhs.getClass() == CountingVariantFilter.class);
        Assert.assertTrue(andFilter.lhs.delegateFilter == snpFilter);
        Assert.assertTrue(andFilter.rhs.getClass() == CountingVariantFilter.class);
        Assert.assertTrue(andFilter.rhs.delegateFilter == biallelicFilter);
    }

}

