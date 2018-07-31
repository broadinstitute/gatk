package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.FisherExactTest;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.Map;


/**
 * Strand bias estimated using Fisher's Exact Test
 *
 * <p>Strand bias is a type of sequencing bias in which one DNA strand is favored over the other, which can result in incorrect evaluation of the amount of evidence observed for one allele vs. the other. The FisherStrand annotation is one of several methods that aims to evaluate whether there is strand bias in the data. It uses Fisher's Exact Test to determine if there is strand bias between forward and reverse strands for the reference or alternate allele.</p>
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
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandBiasBySample.php">StrandBiasBySample</a></b> outputs counts of read depth per allele for each strand orientation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php">StrandOddsRatio</a></b> is an updated form of FisherStrand that uses a symmetric odds ratio calculation.</li>
 * </ul>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Strand bias estimated using Fisher's exact test (FS)")
public final class FisherStrand extends StrandBiasTest implements StandardAnnotation {

    static final double MIN_PVALUE = 1E-320;
    private static final int MIN_COUNT = ARRAY_DIM;
    private static final int MIN_QUAL_FOR_FILTERED_TEST = 17;

    // how large do we want the normalized table to be? (ie, sum of all entries must be smaller that this)
    private static final double TARGET_TABLE_SIZE = 200.0;

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.FISHER_STRAND_KEY);
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        final int[][] tableFromPerSampleAnnotations = getTableFromSamples(genotypes, MIN_COUNT);
        return ( tableFromPerSampleAnnotations != null )? annotationForOneTable(pValueForContingencyTable(tableFromPerSampleAnnotations)) : null;
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromStratifiedContexts(final Map<String, List<PileupElement>> stratifiedContexts,
                                                                            final VariantContext vc){
        final int[][] tableNoFiltering = getPileupContingencyTable(stratifiedContexts, vc.getReference(), vc.getAlternateAlleles(), -1, MIN_COUNT);
        final int[][] tableFiltering = getPileupContingencyTable(stratifiedContexts, vc.getReference(), vc.getAlternateAlleles(), MIN_QUAL_FOR_FILTERED_TEST, MIN_COUNT);
        return annotationForOneTable(Math.max(pValueForContingencyTable(tableFiltering), pValueForContingencyTable(tableNoFiltering)));
    }

    @Override
    protected Map<String, Object> calculateAnnotationFromLikelihoods(final ReadLikelihoods<Allele> likelihoods,
                                                                     final VariantContext vc){
        final int[][] table = getContingencyTable(likelihoods, vc, MIN_COUNT);
        return annotationForOneTable(pValueForContingencyTable(table));
    }

    /**
     * Returns an annotation result given a pValue
     *
     * @param pValue
     * @return a hash map from FS -> phred-scaled pValue
     */
    @VisibleForTesting
    Map<String, Object> annotationForOneTable(final double pValue) {
        return Collections.singletonMap(getKeyNames().get(0), makeValueObjectForAnnotation(pValue));
    }

    public static String makeValueObjectForAnnotation(final int[][] originalTable) {
        return makeValueObjectForAnnotation(pValueForContingencyTable(originalTable));
    }

    public static String makeValueObjectForAnnotation(double pValue) {
        return String.format("%.3f", QualityUtils.phredScaleErrorRate(Math.max(pValue, MIN_PVALUE))); // prevent INFINITYs
    }

    public static Double pValueForContingencyTable(final int[][] originalTable) {
        final int[][] normalizedTable = normalizeContingencyTable(originalTable);
        return FisherExactTest.twoSidedPValue(normalizedTable);
    }

    /**
     * Normalize the table so that the entries are not too large.
     * Note that this method does NOT necessarily make a copy of the table being passed in!
     *
     * @param table  the original table
     * @return a normalized version of the table or the original table if it is already normalized
     */
    private static int[][] normalizeContingencyTable(final int[][] table) {
        final int sum = addExact(table[0][0], table[0][1], table[1][0], table[1][1]);
        if ( sum <= TARGET_TABLE_SIZE * 2 ) {
            return table;
        }

        final double normFactor = sum / TARGET_TABLE_SIZE;

        return new int[][]{
                {(int) (table[0][0] / normFactor), (int) (table[0][1] / normFactor)},
                {(int) (table[1][0] / normFactor), (int) (table[1][1] / normFactor)},
        };
    }

    //Add a bunch of ints, blows up if there's overflow
    private static int addExact(final int... ints){
        int res = ints[0];
        for (int i = 1; i < ints.length; i++) {
            res = Math.addExact(res, ints[i]);
        }
        return res;
    }
}
