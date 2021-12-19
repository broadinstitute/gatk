package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;

import java.util.List;

public class ReblockingGVCFWriter extends GVCFWriter {

    public ReblockingGVCFWriter(final VariantContextWriter underlyingWriter, final List<? extends Number> gqPartitions,
                                final boolean floorBlocks, final CachingIndexedFastaSequenceFile referenceReader,
                                final ReblockingOptions reblockingOptions) {
        super(underlyingWriter, gqPartitions, floorBlocks);
        this.gvcfBlockCombiner = new ReblockingGVCFBlockCombiner(gqPartitions, floorBlocks, referenceReader, reblockingOptions);
    }

    @Override
    public void add(VariantContext vc) {
        gvcfBlockCombiner.submit(vc);
        output();
    }


    /**
     * Return the position of the end of the variant contexts already output by the writer
     * @return may be null
     */
    public Locatable getVcfOutputEnd() {
        final String contig = ((ReblockingGVCFBlockCombiner)gvcfBlockCombiner).getCurrentContig();
        if (contig == null) {
            return null;
        }
        final int position = ((ReblockingGVCFBlockCombiner)gvcfBlockCombiner).getVcfOutputEnd();
        if (position == 0) {
            return null;
        }
        return new SimpleInterval(contig, position, position);
    }


    /**
     *
     * @param vc
     * @return true if the input variant context overlapped the stored (but not output) reference blocks
     */
    public boolean siteOverlapsBuffer(final VariantContext vc) {
        final ReblockingGVCFBlockCombiner combiner = ((ReblockingGVCFBlockCombiner)gvcfBlockCombiner);
        if (combiner.isBufferEmpty()) { return false;}
        return vc.getContig().equals(combiner.getCurrentContig())
            && vc.getStart() <= combiner.getBufferEnd()
            && vc.getStart() >= combiner.getBufferStart(); }  //we need the start too in case there's a tail block from a trimmed deletion, but the deletion still overlaps
}
