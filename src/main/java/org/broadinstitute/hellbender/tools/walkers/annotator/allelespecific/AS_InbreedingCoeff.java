package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.HeterozygosityCalculator;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Allele-specific likelihood-based test for the inbreeding among samples
 *
 * <p>This annotation estimates whether there is evidence of inbreeding in a population. The higher the score, the higher the chance that there is inbreeding.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The calculation is a continuous generalization of the Hardy-Weinberg test for disequilibrium that works well with limited coverage per sample. The output is the F statistic from running the HW test for disequilibrium with PL values. See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of this statistical test.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>The inbreeding coefficient can only be calculated for cohorts containing at least 10 founder samples.</li>
 * <li>This annotation can take a valid pedigree file to specify founders. If not specified, all samples will be considered as founders.</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_InbreedingCoeff.php">InbreedingCoeff</a></b> outputs a version of this annotation that includes all alternate alleles in a single calculation.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_ExcessHet.php">ExcessHet</a></b> estimates excess heterozygosity in a population of samples.</li>
 * </ul>
 *
 */
//TODO: this can't extend InbreedingCoeff because that one is Standard and it would force this to be output all the time; should fix code duplication nonetheless
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Allele-specific likelihood-based test for the consanguinity among samples (AS_InbreedingCoeff)")
public final class AS_InbreedingCoeff extends InfoFieldAnnotation implements AS_StandardAnnotation {

    public static final int MIN_SAMPLES = 10;
    private Set<String> founderIds;    //TODO: either use this or enter a bug report

    public AS_InbreedingCoeff(){
        this(null);
    }

    public AS_InbreedingCoeff(final Set<String> founderIds){
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        this.founderIds = founderIds;
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY); }


    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods ) {
        Utils.nonNull(vc);
        final HeterozygosityCalculator heterozygosityUtils = new HeterozygosityCalculator(vc);

        if (heterozygosityUtils.getSampleCount() < MIN_SAMPLES) {
            return Collections.emptyMap();
        }
        final List<Double> ICvalues = new ArrayList<>();
        for (final Allele a : vc.getAlternateAlleles()) {
            ICvalues.add(calculateIC(vc, a, heterozygosityUtils));
        }
        return Collections.singletonMap(getKeyNames().get(0),  AnnotationUtils.encodeValueList(ICvalues, "%.4f"));
    }

    @VisibleForTesting
    public double calculateIC(final VariantContext vc, final Allele altAllele) {
        //make a new HeterozygosityUtils for each call to reset it
        return calculateIC(vc, altAllele, new HeterozygosityCalculator(vc));
    }

    private double calculateIC(final VariantContext vc, final Allele altAllele, final HeterozygosityCalculator heterozygosityUtils) {
        final int AN = vc.getCalledChrCount();
        final double altAF;

        final double hetCount = heterozygosityUtils.getHetCount(altAllele);

        final double F;
        //shortcut to get a value closer to the non-alleleSpecific value for bialleleics
        if (vc.isBiallelic()) {
            double refAC = heterozygosityUtils.getAlleleCount(vc.getReference());
            double altAC = heterozygosityUtils.getAlleleCount(altAllele);
            double refAF = refAC/(altAC+refAC);
            altAF = 1 - refAF;
            F = 1.0 - (hetCount / (2.0 * refAF * altAF * (double) heterozygosityUtils.getSampleCount())); // inbreeding coefficient
        }
        else {
            //compare number of hets for this allele (and any other second allele) with the expectation based on AFs
            //derive the altAF from the likelihoods to account for any accumulation of fractional counts from non-primary likelihoods,
            //e.g. for a GQ10 variant, the probability of the call will be ~0.9 and the second best call will be ~0.1 so adding up those 0.1s for het counts can dramatically change the AF compared with integer counts
            altAF = heterozygosityUtils.getAlleleCount(altAllele)/ (double) AN;
            F = 1.0 - (hetCount / (2.0 * (1 - altAF) * altAF * (double) heterozygosityUtils.getSampleCount())); // inbreeding coefficient
        }

        return F;
    }
}
