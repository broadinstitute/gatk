package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Annotate with local reference bases.
 *
 * Created by David Benjamin on 3/20/17.
 */
public class ReferenceBases extends InfoFieldAnnotation {
    public static final String REFERENCE_BASES_KEY = "REF_BASES";

    public static final int NUM_BASES_ON_EITHER_SIDE = 10;

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(REFERENCE_BASES_KEY); }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {

        final int basesToDiscardInFront = Math.max(vc.getStart() - ref.getWindow().getStart() - NUM_BASES_ON_EITHER_SIDE, 0);
        final String allBases = new String(ref.getBases());
        final String localBases = allBases.substring(basesToDiscardInFront, basesToDiscardInFront + 2 * NUM_BASES_ON_EITHER_SIDE);
        return Collections.singletonMap(REFERENCE_BASES_KEY, localBases );
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(ReferenceBases.REFERENCE_BASES_KEY, 1, VCFHeaderLineType.String, "local reference bases."));
    }
}