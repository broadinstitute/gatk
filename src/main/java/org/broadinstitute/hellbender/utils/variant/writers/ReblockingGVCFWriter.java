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
     *
     * @return may be null
     */
    public Locatable getVcfOutputEnd() {
        final int position = ((ReblockingGVCFBlockCombiner)gvcfBlockCombiner).getVcfOutputEnd();
        final String contig = ((ReblockingGVCFBlockCombiner)gvcfBlockCombiner).getCurrentContig();
        if (contig == null) {
            return null;
        }
        return new SimpleInterval(contig, position, position);
    }

    public boolean siteOverlapsBuffer(final VariantContext vc) { return vc.getContig().equals(((ReblockingGVCFBlockCombiner)gvcfBlockCombiner).getCurrentContig())
            && vc.getStart() <= ((ReblockingGVCFBlockCombiner)gvcfBlockCombiner).getBufferEnd(); }
}
