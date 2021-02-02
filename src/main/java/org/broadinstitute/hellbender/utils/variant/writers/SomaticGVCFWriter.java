package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.List;

/**
 * Genome-wide VCF writer for somatic (Mutect2) output
 * Merges reference blocks based on TLOD
 */
final public class SomaticGVCFWriter extends GVCFWriter {

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
     * @param lodPartitions     a list of TLOD partitions, this list must be non-empty and every element must be larger than previous element
     */
    public SomaticGVCFWriter(final VariantContextWriter underlyingWriter, final List<Number> lodPartitions) {
        super(underlyingWriter, ImmutableList.of(1, 10, 20), false);
        gvcfBlockCombiner = new SomaticGVCFBlockCombiner(lodPartitions, HomoSapiensConstants.DEFAULT_PLOIDY);
    }

    public SomaticGVCFWriter(final VariantContextWriter underlyingWriter, final List<Number> lodPartitions, final int partitionPrecision) {
        super(underlyingWriter, ImmutableList.of(1, 10, 20), false);
        gvcfBlockCombiner = new SomaticGVCFBlockCombiner(lodPartitions, HomoSapiensConstants.DEFAULT_PLOIDY, partitionPrecision);
    }

    @VisibleForTesting
    protected int convertLODtoInt(final double LOD) {
        return (int)Math.round(LOD * Math.pow(10, ((SomaticGVCFBlockCombiner)gvcfBlockCombiner).partitionPrecision));
    }

}