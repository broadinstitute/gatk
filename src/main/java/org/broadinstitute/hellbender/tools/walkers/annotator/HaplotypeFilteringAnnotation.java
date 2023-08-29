package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Set of annotations meant to be reflective of HaplotypeFiltering operations that were applied in FlowBased HaplotypeCaller.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_FLOW_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_FLOW_ANNOTATORS_SUMMARY,
        summary="Summary of the haplotype filtering steps.")
public class HaplotypeFilteringAnnotation implements JumboInfoAnnotation {


    public HaplotypeFilteringAnnotation() {
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final FeatureContext features,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                        final AlleleLikelihoods<?, Allele> fragmentLikelihoods,
                                        final AlleleLikelihoods<?, Haplotype> haplotypeLikelihoods) {

        final Map<String, Object> result = new HashMap<>();
        result.put(GATKVCFConstants.HAPLOTYPES_BEFORE_FILTERING_KEY, haplotypeLikelihoods.alleles().size());
        result.put(GATKVCFConstants.HAPLOTYPES_FILTERED_KEY, haplotypeLikelihoods.getFilteredHaplotypeCount());

        return result;
    }


    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.HAPLOTYPES_BEFORE_FILTERING_KEY, GATKVCFConstants.HAPLOTYPES_FILTERED_KEY);
    }
}
