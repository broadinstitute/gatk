package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.Utils;

public class SelectVariantsUnitTest extends BaseTest {

    ///////////////////////////////////////////////////////////
    // Tests for maxIndelSize and minIndelSize functionality //
    ///////////////////////////////////////////////////////////

    @DataProvider(name = "MaxMinIndelSize")
    public Object[][] MaxMinIndelSizeTestData() {

        List<Object[]> tests = new ArrayList<>();

        for ( final int size : Arrays.asList(1, 3, 10, 100) ) {
            for ( final int otherSize : Arrays.asList(0, 1) ) {
                for ( final int max : Arrays.asList(0, 1, 5, 50, 100000) ) {
                    for ( final int min : Arrays.asList(0, 1, 5, 50) ) {
                        for (final String op : Arrays.asList("D", "I")) {
                            tests.add(new Object[]{size, otherSize, max, min, op});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MaxMinIndelSize")
    public void maxIndelSizeTest(final int size, final int otherSize, final int max, final int min, final String op) {

        final byte[] largerAllele = Utils.dupBytes((byte) 'A', size+1);
        final byte[] smallerAllele = Utils.dupBytes((byte) 'A', 1);

        final List<Allele> alleles = new ArrayList<>(2);
        final Allele ref = Allele.create(op.equals("I") ? smallerAllele : largerAllele, true);
        final Allele alt = Allele.create(op.equals("D") ? smallerAllele : largerAllele, false);
        alleles.add(ref);
        alleles.add(alt);

        final VariantContext vc = new VariantContextBuilder("test", "1", 10, 10 + ref.length() - 1, alleles).make();

        final boolean hasIndelTooLargeOrSmall = SelectVariants.containsIndelLargerOrSmallerThan(vc, max, min);
        Assert.assertEquals(hasIndelTooLargeOrSmall, size > max || size < min);
    }

    @DataProvider(name="includeMatchingTests")
    public Object[][] getIncludeMatchingTests() {
        return new Object[][] {
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a"), true, new LinkedHashSet<>(Arrays.asList("a")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a"), false, new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("b"), true, Collections.EMPTY_SET },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("b"), false, new LinkedHashSet<>(Arrays.asList("ab", "abc")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "b"), true, new LinkedHashSet<>(Arrays.asList("a")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "b"), false, new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "ab"), true, new LinkedHashSet<>(Arrays.asList("a", "ab")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList("a", "ab"), false, new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")) },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*b.*"), true, Collections.EMPTY_SET },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*b.*"), false, new LinkedHashSet<>(Arrays.asList("ab", "abc") )},
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*"), true, Collections.EMPTY_SET },
                new Object[] { new LinkedHashSet<>(Arrays.asList("a", "ab", "abc")), Arrays.asList(".*"), false, new LinkedHashSet<>(Arrays.asList("a", "ab", "abc") )}
        };
    }

    @Test(dataProvider = "includeMatchingTests")
    public void testIncludeMatching(Set<String> values, Collection<String> filters, boolean exactMatch, Set<String> expected) {
        Set<String> actual = SelectVariants.includeMatching(values, filters, exactMatch);
        Assert.assertEquals(actual, expected);
    }

}