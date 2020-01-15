package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.spark.sql.catalyst.expressions.aggregate.Count;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class CountNsUnitTest {

    @Test
    public void testDoesReadHaveN() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        final String contig = "chr1";

        final int pos = 100;
        final VariantContext vc = new VariantContextBuilder().chr(contig).start(pos).stop(pos).alleles("A","C").make();

        final GATKRead nonOverlappingUpstreamRead = ArtificialReadUtils.createArtificialRead(header, "read", 0, 96, new byte[] {'A','N','G'}, new byte[] {20,20,20}, "3M");
        Assert.assertFalse(CountNs.doesReadHaveN(nonOverlappingUpstreamRead, vc));

        final GATKRead nonOverlappingDownstreamRead = ArtificialReadUtils.createArtificialRead(header, "read", 0, 101, new byte[] {'A','N','G'}, new byte[] {20,20,20}, "3M");
        Assert.assertFalse(CountNs.doesReadHaveN(nonOverlappingDownstreamRead, vc));

        final GATKRead spanningDeletionRead = ArtificialReadUtils.createArtificialRead(header, "read", 0, 95, new byte[] {'N','N','N','N','N','N'}, new byte[] {20,20,20,20,20,20}, "3M10D3M");
        Assert.assertFalse(CountNs.doesReadHaveN(spanningDeletionRead, vc));

        final GATKRead notN = ArtificialReadUtils.createArtificialRead(header, "read", 0, 99, new byte[] {'A','C','G'}, new byte[] {20,20,20}, "3M");
        Assert.assertFalse(CountNs.doesReadHaveN(notN, vc));

        final GATKRead yesN = ArtificialReadUtils.createArtificialRead(header, "read", 0, 99, new byte[] {'A','N','G'}, new byte[] {20,20,20}, "3M");
        Assert.assertTrue(CountNs.doesReadHaveN(yesN, vc));

    }
}