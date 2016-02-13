package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
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
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        Utils.nonNull(vc);
        final GenotypesContext genotypes = (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(founderIds);
        if (genotypes == null || genotypes.size() < MIN_SAMPLES || !vc.isVariant()) {
            return null;
        }
        final Pair<Integer, Double> sampleCountCoeff = calculateIC(vc, genotypes);
        final int sampleCount = sampleCountCoeff.getLeft();
        final double F = sampleCountCoeff.getRight();
        if (sampleCount < MIN_SAMPLES) {
            logger.warn("Annotation will not be calculated, must provide at least " + MIN_SAMPLES + " samples");
            return null;
        }
        return Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", F));
    }

    @VisibleForTesting
    Pair<Integer, Double> calculateIC(final VariantContext vc, final GenotypesContext genotypes) {

        final boolean doMultiallelicMapping = !vc.isBiallelic();

        int idxAA = 0, idxAB = 1, idxBB = 2;

        double refCount = 0.0;
        double hetCount = 0.0;
        double homCount = 0.0;
        int sampleCount = 0; // number of samples that have likelihoods

        for ( final Genotype g : genotypes ) {
            if ( g.isCalled() && g.hasLikelihoods() && g.getPloidy() == 2){  // only work for diploid samples
                sampleCount++;
            } else {
                continue;
            }
            final double[] normalizedLikelihoods = MathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
            if (doMultiallelicMapping){
                if (g.isHetNonRef()) {
                    //all likelihoods go to homCount
                    homCount++;
                    continue;
                }

                //get alternate allele for each sample
                final Allele a1 = g.getAllele(0);
                final Allele a2 = g.getAllele(1);
                if (a2.isNonReference()) {
                    final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a2);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
                //I expect hets to be reference first, but there are no guarantees (e.g. phasing)
                else if (a1.isNonReference()) {
                    final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a1);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
            }

            refCount += normalizedLikelihoods[idxAA];
            hetCount += normalizedLikelihoods[idxAB];
            homCount += normalizedLikelihoods[idxBB];
        }

        /**
         * Note: all that likelihood normalization etc may have accumulated some error.
         * We smooth it out my roudning the numbers to integers before the final computation.
         */
        refCount = Math.round(refCount);
        hetCount = Math.round(hetCount);
        homCount = Math.round(homCount);

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