package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class DownsampledFlag extends InfoFieldAnnotation implements StandardHCAnnotation {
    /**
     * Computes the annotation for the given variant and the likelihoods per read.
     * Returns a map from annotation keys to values (may be empty if no annotation is to be added).
     *
     * @param ref         Reference context, may be null
     * @param vc          Variant to be annotated. Not null.
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     */
    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        if (likelihoods.getCoverageDownsamplingHistory()) {
            return Collections.singletonMap(getKeyNames().get(0), null);
        } else {
            return null;
        }
    }

    /**
     * Return the keys
     */
    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.DOWNSAMPLED_KEY);
    }
}
