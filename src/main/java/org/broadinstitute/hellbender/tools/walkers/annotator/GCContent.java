package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

/**
 * Created by tsato on 6/30/17.
 */
public class GCContent extends InfoFieldAnnotation {
    private static final String KEY = "GC";
    private static final int k = 100;

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        // TODO: is this right?
        return Collections.singletonList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, ReadLikelihoods<Allele> likelihoods) {
        SimpleInterval oldWindow = ref.getWindow();
        // TODO: make interval precise
        ref.setWindow(k/2, k/2);
        final int C = 67;
        final int G = 71;
        final byte[] bases = ref.getBases();
        final long gccount = IntStream.range(0, bases.length).map(i -> bases[i]).filter(j -> j == C || j == G).count();
        final double gcfraction = (double) gccount / bases.length;
        // TODO: restore window?
        return Collections.singletonMap(getKeyNames().get(0), String.format("%.4f", gcfraction));
    }
}
