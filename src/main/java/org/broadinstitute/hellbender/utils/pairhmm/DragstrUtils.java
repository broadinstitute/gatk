package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.aeonbits.owner.util.Collections;

import java.util.Collection;

public class DragstrUtils {

    public static DragstrReadSTRAnalizer repeatPeriodAndCounts(final int maxSequenceLength, final int maxPeriod, final boolean considerUpstream) {
        return new DragstrReadSTRAnalizer(maxSequenceLength, maxPeriod, considerUpstream);
    }

    public static DragstrReadSTRAnalizer repeatPeriodAndCounts(final byte[] sequence, final int maxPeriod, final boolean considerUpstream) {
        final DragstrReadSTRAnalizer result = new DragstrReadSTRAnalizer(sequence.length, maxPeriod, considerUpstream);
        result.load(sequence);
        return result;
    }

    public static DragstrReadSTRAnalizer repeatPeriodAndCounts(final byte[] sequence, final int start, final int stop, final int maxPeriod, final boolean considerUpstream) {
        final DragstrReadSTRAnalizer result = new DragstrReadSTRAnalizer(sequence.length, maxPeriod, considerUpstream);
        result.load(sequence, start, stop);
        return result;
    }

    public static Collection<? extends VCFHeaderLine> vcfHeaderLines() {
        return Collections.list(new VCFInfoHeaderLine(DragstrConstants.DRAGSTRINFO_KEY, 2, VCFHeaderLineType.Integer, "Indicates the period and repeat count"),
                                new VCFInfoHeaderLine(DragstrConstants.DRAGSTRPARAMS_KEY, 3, VCFHeaderLineType.Float, "Parameeters used (GOP, GCP, API)"));

    }

    public static VariantContext annotate(VariantContext annotatedCall, final DragstrParams dragstrParams, final int period, final int repeats) {
        final VariantContextBuilder builder = new VariantContextBuilder(annotatedCall);
        final double gop = dragstrParams.gop(period, repeats);
        final double gcp = dragstrParams.gcp(period, repeats);
        final double api = dragstrParams.api(period, repeats);
        builder.attribute(DragstrConstants.DRAGSTRINFO_KEY, new int[] {period, repeats});
        builder.attribute(DragstrConstants.DRAGSTRPARAMS_KEY, new String[] {String.format("%.1f", gop),String.format("%.1f", gcp), String.format("%.1f", api)});
        return builder.make();
    }

}
