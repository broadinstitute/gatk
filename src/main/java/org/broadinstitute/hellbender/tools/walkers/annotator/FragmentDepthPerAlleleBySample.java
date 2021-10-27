package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Fragment depth of coverage of each allele per sample
 *
 * <p>This annotation is identical to DepthPerAlleleBySample except that allele support is considered per read pair, not per individual read.</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Depth of coverage of each allele per sample (AD)")
public final class FragmentDepthPerAlleleBySample implements JumboGenotypeAnnotation, StandardMutectAnnotation {

    @Override
    public void annotate(final ReferenceContext ref,
                         final FeatureContext features,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> readLikelihoods,
                         final AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods,
                         final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods) {
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        if ( g == null || !g.isCalled() || fragmentLikelihoods == null) {
            return;
        }
        final Set<Allele> alleles = new LinkedHashSet<>(vc.getAlleles());

        // make sure that there's a meaningful relationship between the alleles in the likelihoods and our VariantContext
        Utils.validateArg(fragmentLikelihoods.alleles().containsAll(alleles), () -> "VC alleles " + alleles + " not a  subset of AlleleLikelihoods alleles " + fragmentLikelihoods.alleles());


        final int[] alleleDepths = DepthPerAlleleBySample.annotateWithLikelihoods(vc, g, alleles, fragmentLikelihoods);
        gb.attribute(GATKVCFConstants.FRAGMENT_ALLELE_DEPTHS, alleleDepths);
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.FRAGMENT_ALLELE_DEPTHS); }
}
