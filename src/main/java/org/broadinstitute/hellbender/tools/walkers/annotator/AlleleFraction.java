package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Created by gauthier on 8/15/18.
 * Bug Laura to fill this in
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Variant allele fraction for a genotype")

public final class AlleleFraction extends GenotypeAnnotation implements StandardMutectAnnotation {
    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( g == null || !g.isCalled() || g.hasExtendedAttribute(getKeyNames().get(0))) {  //don't overwrite AF based on Bayesian estimate if it already exists
            return;
        }

        for ( final Genotype genotype : genotypes ) {
           if ( genotype.hasAD() ) {
                final int[] AD = genotype.getAD();
                final double[] allAlleleFractions = MathUtils.normalizeFromRealSpace(Arrays.stream(AD).mapToDouble(x -> x*1.0).toArray());
                gb.attribute(getKeyNames().get(0), Arrays.copyOfRange(allAlleleFractions, 1, allAlleleFractions.length)); //omit the first entry of the array corresponding to the reference
            }
            // if there is no AD value calculate it now using likelihoods
            else if (likelihoods != null) {
                DepthPerAlleleBySample adCalc = new DepthPerAlleleBySample();
                final int[] AD = adCalc.annotateWithLikelihoods(vc, g, new LinkedHashSet<>(vc.getAlleles()), likelihoods);
                final double[] allAlleleFractions = MathUtils.normalizeFromRealSpace(Arrays.stream(AD).mapToDouble(x -> x*1.0).toArray());
                gb.attribute(getKeyNames().get(0), Arrays.copyOfRange(allAlleleFractions, 1, allAlleleFractions.length)); //omit the first entry of the array corresponding to the reference
            }
        }
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.ALLELE_FRACTION_KEY);}

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getFormatLine(getKeyNames().get(0)));
    }
}
