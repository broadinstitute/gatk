package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * Likelihood-based test for the inbreeding among samples
 *
 * <p>This annotation estimates whether there is evidence of inbreeding in a population. The higher the score, the higher the chance that there is inbreeding.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The calculation is a continuous generalization of the Hardy-Weinberg test for disequilibrium that works well with limited coverage per sample. The output is a Phred-scaled p-value derived from running the HW test for disequilibrium with PL values. See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of this statistical test.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>The Inbreeding Coefficient annotation can only be calculated for cohorts containing at least 10 founder samples.</li>
 * <li>The Inbreeding Coefficient annotation can only be calculated for diploid samples.</li>
 * <li>This annotation is used in variant recalibration, but may not be appropriate for that purpose if the cohort being analyzed contains many closely related individuals.</li>
 * </ul>
 *
 */
public final class InbreedingCoeff extends InfoFieldAnnotation implements StandardAnnotation {

    private static final Logger logger = LogManager.getLogger(InbreedingCoeff.class);
    private static final int MIN_SAMPLES = 10;
    private final Set<String> founderIds;

    public InbreedingCoeff(){
        this(null);
    }

    public InbreedingCoeff(final Set<String> founderIds){
        //If available, get the founder IDs and cache them. the IC will only be computed on founders then.
        this.founderIds = founderIds;
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        final GenotypesContext genotypes = (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(founderIds);
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
        return Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", F));
    }

    @VisibleForTesting
    Pair<Integer, Double> calculateIC(final VariantContext vc, final GenotypesContext genotypes) {
        final GenotypeCounts t = GenotypeUtils.computeDiploidGenotypeCounts(vc, genotypes);

        final int refCount = t.getRefs();
        final int hetCount = t.getHets();
        final int homCount = t.getHoms();
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

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() { return Collections.singletonList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0))); }
}