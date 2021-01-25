package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.Map;

public class AS_VarDP extends InfoFieldAnnotation implements ReducibleAnnotation, AlleleSpecificAnnotation {
    @Override
    public String getPrimaryRawKey() {
        return null;
    }

    @Override
    public boolean hasSecondaryRawKeys() {
        return false;
    }

    @Override
    public List<String> getSecondaryRawKeys() {
        return null;
    }

    @Override
    public Map<String, Object> annotateRawData(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        return null;
    }

    @Override
    public Map<String, Object> combineRawData(List<Allele> allelesList, List<ReducibleAnnotationData<?>> listOfRawData) {
        return null;
    }

    @Override
    public Map<String, Object> finalizeRawData(VariantContext vc, VariantContext originalVC) {
        return null;
    }

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        return null;
    }

    @Override
    public List<String> getKeyNames() {
        return null;
    }
}
