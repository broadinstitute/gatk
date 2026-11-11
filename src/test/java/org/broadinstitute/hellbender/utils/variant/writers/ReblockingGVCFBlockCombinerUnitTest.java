package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class ReblockingGVCFBlockCombinerUnitTest extends GATKBaseTest {
    static final String CONTIG = "chr20";
    private static final Allele REF = Allele.create("A", true);

    // No reference: these tests never move the start of a block, which is all the combiner reads it for
    private static ReblockingGVCFBlockCombiner newCombiner() {
        return new ReblockingGVCFBlockCombiner(Arrays.asList(20, 40), false, null, new ReblockingOptions());
    }

    static VariantContext refBlock(final int start, final int end, final int gq) {
        final Genotype g = new GenotypeBuilder("sample", Arrays.asList(REF, REF)).GQ(gq).DP(10).PL(new int[]{0, gq, 10 * gq}).make();
        return new VariantContextBuilder("test", CONTIG, start, end, Arrays.asList(REF, Allele.NON_REF_ALLELE))
                .attribute(VCFConstants.END_KEY, end).genotypes(g).make();
    }

    @Test
    public void testFlushEmitsTheBlocksEndingBeforeThePosition() {
        final ReblockingGVCFBlockCombiner combiner = newCombiner();
        combiner.submit(refBlock(100, 199, 50));
        combiner.submit(refBlock(200, 299, 10));
        combiner.submit(refBlock(300, 399, 50));
        Assert.assertFalse(combiner.hasFinalizedItems());

        combiner.flushBlocksEndingBefore(CONTIG, 350);
        Assert.assertEquals(combiner.getBufferStart(), 300);
        Assert.assertEquals(combiner.getVcfOutputEnd(), 299);
        // The band combiner finishes the first block when the second block, in another band, reaches it
        final List<VariantContext> flushed = combiner.consumeFinalizedItems();
        Assert.assertEquals(flushed.size(), 1);
        Assert.assertEquals(flushed.get(0).getStart(), 100);
        Assert.assertEquals(flushed.get(0).getEnd(), 199);

        combiner.signalEndOfInput();
        final List<VariantContext> rest = combiner.consumeFinalizedItems();
        Assert.assertEquals(rest.size(), 2);
        Assert.assertEquals(rest.get(0).getStart(), 200);
        Assert.assertEquals(rest.get(1).getStart(), 300);
    }

    @DataProvider
    public Object[][] flushPositions() {
        // the position to flush before, and whether the block 100-199 stays buffered
        return new Object[][]{{198, true}, {199, true}, {200, false}, {201, false}};
    }

    @Test(dataProvider = "flushPositions")
    public void testFlushKeepsBlocksEndingAtOrAfterThePosition(final int position, final boolean staysBuffered) {
        final ReblockingGVCFBlockCombiner combiner = newCombiner();
        combiner.submit(refBlock(100, 199, 50));

        combiner.flushBlocksEndingBefore(CONTIG, position);
        Assert.assertEquals(combiner.isBufferEmpty(), !staysBuffered);
        Assert.assertEquals(combiner.getVcfOutputEnd(), staysBuffered ? 0 : 199);
    }

    // Emptying the buffer resets its end, as writing the whole buffer does
    @Test
    public void testFlushOfTheWholeBufferResetsItsEnd() {
        final ReblockingGVCFBlockCombiner combiner = newCombiner();
        combiner.submit(refBlock(100, 199, 50));
        Assert.assertEquals(combiner.getBufferEnd(), 199);

        combiner.flushBlocksEndingBefore(CONTIG, 300);
        Assert.assertTrue(combiner.isBufferEmpty());
        Assert.assertEquals(combiner.getBufferEnd(), 0);
    }

    @Test
    public void testFlushOnAnotherContigDoesNothing() {
        final ReblockingGVCFBlockCombiner combiner = newCombiner();
        combiner.flushBlocksEndingBefore(CONTIG, 1000);  // nothing submitted yet
        combiner.submit(refBlock(100, 199, 50));

        combiner.flushBlocksEndingBefore("chr21", 1000);
        Assert.assertEquals(combiner.getBufferStart(), 100);
        Assert.assertEquals(combiner.getVcfOutputEnd(), 0);
    }

    // submit() keeps every buffered block after the output, so this sets the output end directly to break that
    @Test(expectedExceptions = IllegalStateException.class)
    public void testFlushRejectsABlockEndingBeforeTheOutput() {
        final ReblockingGVCFBlockCombiner combiner = newCombiner();
        combiner.submit(refBlock(100, 199, 50));
        combiner.setVcfOutputEnd(250);

        combiner.flushBlocksEndingBefore(CONTIG, 300);
    }
}
