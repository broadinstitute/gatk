package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.aeonbits.owner.util.Collections;

import java.util.Collection;

public class DragstrUtils {

    /**
     * Returns VCF header lines required for annotations supplied by {@link #annotateVariantContextWithDragstrParametersUsed}.
     * @return never {@code null}.
     */
    public static Collection<? extends VCFHeaderLine> vcfHeaderLines() {
        return Collections.list(new VCFInfoHeaderLine(DragstrConstants.DRAGSTRINFO_KEY, 2, VCFHeaderLineType.Integer, "Indicates the period and repeat count"),
                                new VCFInfoHeaderLine(DragstrConstants.DRAGSTRPARAMS_KEY, 3, VCFHeaderLineType.Float, "Parameeters used (GOP, GCP, API)"));

    }

    /**
     * Annotates a variant context with some information on the STR at that site.
     * @param vc the target STR
     * @param dragstrParams DRAGstr model parameters collection.
     * @param period the STR detected period at that position.
     * @param repeats the STR repeats length at that position.
     * @return the updated variant context.
     */
    public static VariantContext annotateVariantContextWithDragstrParametersUsed(final VariantContext vc, final DragstrParams dragstrParams, final int period, final int repeats) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final double gop = dragstrParams.gop(period, repeats);
        final double gcp = dragstrParams.gcp(period, repeats);
        final double api = dragstrParams.api(period, repeats);
        builder.attribute(DragstrConstants.DRAGSTRINFO_KEY, new int[] {period, repeats});
        builder.attribute(DragstrConstants.DRAGSTRPARAMS_KEY, new String[] {String.format("%.1f", gop),String.format("%.1f", gcp), String.format("%.1f", api)});
        return builder.make();
    }

}
