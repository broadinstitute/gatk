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
 * Variant allele fraction for each sample.
 *
 * <p>This annotation describes the proportion of the sample's reads that support the variant allele(s). It uses only reads that are actually considered informative by HaplotypeCaller (HC) or Mutect2, using pre-read likelihoods that are produced internally by HC/Mutect2.</p>
 * <p>In this context, an informative read is defined as one that allows the allele it carries to be easily distinguished. In contrast, a read might be considered uninformative if, for example, it only partially overlaps a short tandem repeat and it is not clear whether the read contains the reference allele or an extra repeat.</p>
 *
 * <p>See the method documentation on <a href="http://www.broadinstitute.org/gatk/guide/article?id=4721">using coverage information</a> for important interpretation details.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>If a genotype is already annotated with allele fraction (as by the SomaticGenotypingEngine if `--get-af-from-ad` is not specified), the value will not be recalculated.</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php">DepthPerAlleleBySample</a></b> calculates depth of coverage for each allele per sample (AD).</li>
 * </ul>
 *
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
                final double[] allAlleleFractions = MathUtils.normalizeFromRealSpace(Arrays.stream(AD).mapToDouble(x -> x).toArray());
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
