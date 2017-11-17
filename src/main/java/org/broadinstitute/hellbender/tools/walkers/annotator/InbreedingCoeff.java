package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;


/**
 * Likelihood-based test for the consanguinuity among samples
 *
 * <p>This annotation estimates whether there is evidence of consanguinuity in a population. The higher the score, the
 * higher the chance that some samples are related. If samples are known to be related, a pedigree file can be provided so
 * that the calculation is only performed on founders and offspring are excluded.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The calculation is a continuous generalization of the Hardy-Weinberg test for disequilibrium that works well with limited coverage per sample. The output is a Phred-scaled p-value derived from running the HW test for disequilibrium with PL values. See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of this statistical test.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>The Inbreeding Coefficient annotation can only be calculated for cohorts containing at least 10 founder samples.</li>
 * <li>The Inbreeding Coefficient annotation can only be calculated for diploid samples.</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <p><b>ExcessHet</b> also describes the heterozygosity of the called samples, giving a probability of excess heterozygosity being observed</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Likelihood-based test for the consanguinity among samples (InbreedingCoeff)")
public final class InbreedingCoeff extends PedigreeAnnotation implements StandardAnnotation {

    private static final Logger logger = LogManager.getLogger(InbreedingCoeff.class);
    private static final int MIN_SAMPLES = 10;
    private static final boolean ROUND_GENOTYPE_COUNTS = false;

    public InbreedingCoeff(){
        super((Set<String>) null);
    }

    public InbreedingCoeff(final Set<String> founderIds){
        super(founderIds);
     }

    public InbreedingCoeff(final File pedigreeFile){
        super(pedigreeFile);
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        final GenotypesContext genotypes = getFounderGenotypes(vc);
        if (genotypes == null || genotypes.size() < MIN_SAMPLES || !vc.isVariant()) {
            return Collections.emptyMap();
        }
        final Pair<Integer, Double> sampleCountCoeff = calculateIC(vc, genotypes);
        final int sampleCount = sampleCountCoeff.getLeft();
        final double F = sampleCountCoeff.getRight();
        if (sampleCount < MIN_SAMPLES) {
            logger.warn("Annotation will not be calculated, must provide at least " + MIN_SAMPLES + " samples");
            return Collections.emptyMap();
        }
        return Collections.singletonMap(getKeyNames().get(0), String.format("%.4f", F));
    }

    @VisibleForTesting
    static Pair<Integer, Double> calculateIC(final VariantContext vc, final GenotypesContext genotypes) {
        final GenotypeCounts t = GenotypeUtils.computeDiploidGenotypeCounts(vc, genotypes, ROUND_GENOTYPE_COUNTS);

        final double refCount = t.getRefs();
        final double hetCount = t.getHets();
        final double homCount = t.getHoms();
        // number of samples that have likelihoods
        final int sampleCount = (int) genotypes.stream().filter(g-> GenotypeUtils.isDiploidWithLikelihoods(g)).count();

        final double p = ( 2.0 * refCount + hetCount ) / ( 2.0 * (refCount + hetCount + homCount) ); // expected reference allele frequency
        final double q = 1.0 - p; // expected alternative allele frequency
        final double expectedHets = 2.0 * p * q * sampleCount; //numbers of hets that would be expected based on the allele frequency (asuming Hardy Weinberg Equilibrium)
        final double F = 1.0 - ( hetCount / expectedHets ); // inbreeding coefficient

        return Pair.of(sampleCount, F);
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY); }
}