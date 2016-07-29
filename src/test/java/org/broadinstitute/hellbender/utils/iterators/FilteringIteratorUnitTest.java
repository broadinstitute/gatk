package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.function.Predicate;

public class FilteringIteratorUnitTest extends BaseTest {

    private Iterator<GATKRead> makeReadsIterator() {
        return Arrays.asList(
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("101M")),
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("50M")),
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("101M")),
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("30M")),
            ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("101M"))
        ).iterator();
    }

    private Iterator<VariantContext> makeVariantIterator() {
        final VariantContextBuilder builder = new VariantContextBuilder().source("test").chr("1").alleles(
            Arrays.asList(Allele.create("A", true), Allele.create("T")));
        int pos = 1;
        return Arrays.asList(
            // 3 biallelic
            builder.start(pos).stop(pos++).make(),
            builder.start(pos).stop(pos++).make(),
            builder.start(pos).stop(pos++).make(),
            // 2 no biallelic
            builder.start(pos).stop(pos++).alleles(Collections.singleton(Allele.create("C", true))).make(),
            builder.start(pos).stop(pos).make()
        ).iterator();
    }

    @DataProvider(name = "FilteringIteratorTestData")
    public Object[][] filteringIteratorTestData() {
        // read filters
        final ReadFilter allowNoReadsFilter = read -> false;
        final ReadFilter allowLongReadsFilter = read -> read.getLength() > 100;
        // variant filters
        final VariantFilter allowNoVariantFilter = variant -> false;
        final VariantFilter allowNoBiallelic = VariantContext::isBiallelic;

        return new Object[][] {
            // testing read filters
            {makeReadsIterator(), ReadFilterLibrary.ALLOW_ALL_READS, 5},
            {makeReadsIterator(), allowNoReadsFilter, 0},
            {makeReadsIterator(), allowLongReadsFilter, 3},
            {makeVariantIterator(), VariantFilterLibrary.ALLOW_ALL_VARIANTS, 5 },
            {makeVariantIterator(), allowNoVariantFilter, 0 },
            {makeVariantIterator(), allowNoBiallelic, 3 }};
    }

    @Test(dataProvider = "FilteringIteratorTestData")
    public void testFilteringIterator(final Iterator<Object> iter, final Predicate<Object> filter, final int expectedReadCount) {
        final FilteringIterator<Object> filteringIterator = new FilteringIterator<>(iter, filter);
        int count = 0;
        for (final Object obj : filteringIterator) {
            ++count;
        }
        Assert.assertEquals(count, expectedReadCount, "Wrong number of reads from the ReadFilteringIterator");
    }
}
