package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

import static org.broadinstitute.hellbender.utils.variant.writers.ReblockingGVCFBlockCombinerUnitTest.CONTIG;
import static org.broadinstitute.hellbender.utils.variant.writers.ReblockingGVCFBlockCombinerUnitTest.refBlock;

public class ReblockingGVCFWriterUnitTest extends GATKBaseTest {

    // The flush writes the blocks it finishes to the underlying writer, without waiting for the next record
    @Test
    public void testFlushWritesTheFinishedBlocks() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final ReblockingGVCFWriter writer = new ReblockingGVCFWriter(mockWriter, Arrays.asList(20, 40), false, null, new ReblockingOptions());
        writer.add(refBlock(100, 199, 50));
        writer.add(refBlock(200, 299, 10));
        writer.add(refBlock(300, 399, 50));
        Assert.assertTrue(mockWriter.emitted.isEmpty());

        writer.flushBlocksEndingBefore(new SimpleInterval(CONTIG, 350, 350));
        // The band combiner finishes the first block when the second block, in another band, reaches it
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        final VariantContext written = mockWriter.emitted.get(0);
        Assert.assertEquals(written.getStart(), 100);
        Assert.assertEquals(written.getEnd(), 199);
    }
}
