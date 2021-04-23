package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.aeonbits.owner.util.Collections;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;

import java.util.Collection;

/**
 * Provides support for Dragstr related annotation on variant-contexts.
 */
final class DragstrVariantContextAnnotations {

    public static final String DRAGSTRINFO_KEY = "DRAGstrInfo";
    public static final String DRAGSTRPARAMS_KEY = "DRAGstrParams";

    /**
     * Returns VCF header lines required for annotations supplied by {@link #annotateVariantContextWithDragstrParametersUsed}.
     * @return never {@code null}.
     */
    static Collection<? extends VCFHeaderLine> vcfHeaderLines() {
        return Collections.list(new VCFInfoHeaderLine(DRAGSTRINFO_KEY, 2, VCFHeaderLineType.Integer, "Indicates the period and repeat count"),
                                new VCFInfoHeaderLine(DRAGSTRPARAMS_KEY, 3, VCFHeaderLineType.Float, "Parameters used (GOP, GCP, API)"));
    }

    /**
     * Annotates a variant context with some information on the STR at that site.
     * @param vc the target STR
     * @param dragstrParams DRAGstr model parameters collection.
     * @param period the STR detected period at that position.
     * @param repeats the STR repeats length at that position.
     * @return the updated variant context.
     */
    static VariantContext annotateVariantContextWithDragstrParametersUsed(final VariantContext vc, final DragstrParams dragstrParams, final int period, final int repeats) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final double gop = dragstrParams.gop(period, repeats);
        final double gcp = dragstrParams.gcp(period, repeats);
        final double api = dragstrParams.api(period, repeats);
        builder.attribute(DRAGSTRINFO_KEY, new int[] {period, repeats});
        builder.attribute(DRAGSTRPARAMS_KEY, new String[] {String.format("%.1f", gop),String.format("%.1f", gcp), String.format("%.1f", api)});
        return builder.make();
    }
}
