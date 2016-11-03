package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * Genome-wide VCF writer
 */
public final class GVCFWriter implements VariantContextWriter {

    public final static String GVCF_BLOCK = "GVCFBlock";

    /** Where we'll ultimately write our VCF records */
    private final VariantContextWriter underlyingWriter;
    private final GVCFBlockCombiner gvcfBlockCombiner;

    /**
     * Create a new GVCF writer
     *
     * Should be a non-empty list of boundaries.  For example, suppose this variable is
     *
     * [A, B, C]
     *
     * We would partition our hom-ref sites into the following bands:
     *
     * X < A
     * A <= X < B
     * B <= X < C
     * X >= C
     *
     * @param underlyingWriter the ultimate destination of the GVCF records
     * @param gqPartitions     a list of GQ partitions, this list must be non-empty and every element must be larger than previous element
     * @param defaultPloidy    the assumed ploidy for input variant context without one.
     */
    public GVCFWriter(final VariantContextWriter underlyingWriter, final List<Integer> gqPartitions, final int defaultPloidy) {
        this(underlyingWriter, new GVCFBlockCombiner(gqPartitions, defaultPloidy));
    }

    public GVCFWriter(final VariantContextWriter underlyingWriter, GVCFBlockCombiner combiner) {
        this.underlyingWriter = Utils.nonNull(underlyingWriter);
        this.gvcfBlockCombiner = Utils.nonNull(combiner);
    }

    /**
     * Write the VCF header
     *
     * Adds standard GVCF fields to the header
     *
     * @param header a non-null header
     */
    @Override
    public void writeHeader(VCFHeader header) {
        gvcfBlockCombiner.addRangesToHeader(header);
        underlyingWriter.writeHeader(header);
    }

    /**
     * Close this GVCF writer.  Finalizes any pending hom-ref blocks and emits those to the underlyingWriter as well
     */
    @Override
    public void close() {
        try {
            gvcfBlockCombiner.signalEndOfInput();
            output();
        } finally {
            underlyingWriter.close();
        }
    }

    @Override
    public boolean checkError() {
        return underlyingWriter.checkError();
    }

    /**
     * Add a VariantContext to this writer for emission
     *
     * Requires that the VC have exactly one genotype
     *
     * @param vc a non-null VariantContext
     */
    @Override
    public void add(VariantContext vc) {
        gvcfBlockCombiner.submit(vc);
        output();
    }

    private void output() {
        if (gvcfBlockCombiner.hasFinalizedItems()) {
            gvcfBlockCombiner.consumeFinalizedItems().forEach(underlyingWriter::add);
        }
    }

    @Override
    public void setHeader(VCFHeader header) {
        gvcfBlockCombiner.addRangesToHeader(header);
        underlyingWriter.setHeader(header);
    }
}
