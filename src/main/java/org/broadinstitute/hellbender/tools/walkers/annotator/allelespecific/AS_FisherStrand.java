package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.tools.walkers.annotator.FisherStrand;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasTest;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Allele-specific strand bias estimated using Fisher's Exact Test
 *
 * * <p>Strand bias is a type of sequencing bias in which one DNA strand is favored over the other, which can result in incorrect evaluation of the amount of evidence observed for one allele vs. the other.</p>
 *
 * <p>The AS_FisherStrand annotation is one of several methods that aims to evaluate whether there is strand bias in the data. It uses Fisher's Exact Test to determine if there is strand bias between forward and reverse strands for the reference or alternate allele, and does so separately for each alternate allele.</p>
 * <p>The output is a Phred-scaled p-value. The higher the output value, the more likely there is to be bias. More bias is indicative of false positive calls.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of this application of Fisher's Exact Test.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>The FisherStrand test may not be calculated for certain complex indel cases or for multi-allelic sites.</li>
 *     <li>FisherStrand is best suited for low coverage situations. For testing strand bias in higher coverage situations, see the StrandOddsRatio annotation.</li>
 * </ul>
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_AS_FisherStrand.php">AS_FisherStrand</a></b> outputs a version of this annotation that includes all alternate alleles in a single calculation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandBiasBySample.php">StrandBiasBySample</a></b> outputs counts of read depth per allele for each strand orientation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php">StrandOddsRatio</a></b> is an updated form of FisherStrand that uses a symmetric odds ratio calculation.</li>
 * </ul>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Allele-specific strand bias estimated using Fisher's exact test (AS_FS)")
public class AS_FisherStrand extends AS_StrandBiasTest implements AS_StandardAnnotation {

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.AS_FISHER_STRAND_KEY);
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoods(final ReadLikelihoods<Allele> likelihoods,
                                                                     final VariantContext vc) {
        // either SNP with no alignment context, or indels: per-read likelihood map needed
        final int[][] table = StrandBiasTest.getContingencyTable(likelihoods, vc, MIN_COUNT);
        return table == null ? null : annotationForOneTable(FisherStrand.pValueForContingencyTable(table));
    }

    /**
     * Returns an annotation result given a pValue
     *
     * @return a hash map from FS -> phred-scaled pValue
     */
    private Map<String, Object> annotationForOneTable(final double pValue) {
        return Collections.singletonMap(getKeyNames().get(0), FisherStrand.makeValueObjectForAnnotation(pValue));
    }

    @Override
    protected Map<Allele,Double> calculateReducedData(AlleleSpecificAnnotationData<List<Integer>> combinedData) {
        final Map<Allele,Double> annotationMap = new HashMap<>();
        final Map<Allele,List<Integer>> perAlleleData = combinedData.getAttributeMap();
        final List<Integer> refStrandCounts = perAlleleData.get(combinedData.getRefAllele());
        for (final Allele a : perAlleleData.keySet()) {
            if(!a.equals(combinedData.getRefAllele(),true)) {
                final List<Integer> altStrandCounts = combinedData.getAttribute(a);
                final int[][] refAltTable = new int[][]{new int[]{refStrandCounts.get(0), refStrandCounts.get(1)}, new int[]{altStrandCounts.get(0), altStrandCounts.get(1)}};
                annotationMap.put(a, QualityUtils.phredScaleErrorRate(Math.max(FisherStrand.pValueForContingencyTable(refAltTable), MIN_PVALUE)));
            }
        }
        return annotationMap;
    }

}
