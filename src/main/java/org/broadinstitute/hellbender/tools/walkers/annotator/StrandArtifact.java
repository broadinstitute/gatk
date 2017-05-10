package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.annotator.StrandArtifact.StrandArtifactZ.*;

/**
 * Annotations for strand artifact filter
 *
 * Created by tsato on 4/19/17.
 */
public class StrandArtifact extends GenotypeAnnotation implements StandardSomaticAnnotation {
    public static final String POSTERIOR_PROBABILITIES_KEY = "SA_POST_PROB";
    public static final String MAP_ALLELE_FRACTIONS_KEY = "SA_MAP_AF";

    // pseudocounts for the beta distribution over epsilon
    // alpha > 0 and beta > 0. alpha = beta = 1 gives us the flat prior
    // give more prior weight to beta (i.e. pseudocount for tails) so that the peak shifts towards 0
    private static final int ALPHA = 1;
    private static final int BETA = 6;
    private static final BetaDistribution betaPrior = new BetaDistribution(null, ALPHA, BETA);

    // the latent variable z in the strand artifact filter model
    public enum StrandArtifactZ {
        ART_FWD, ART_REV, NO_ARTIFACT
    }

    private static final EnumMap<StrandArtifactZ, Double> pi = new EnumMap<>(StrandArtifactZ.class);

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(POSTERIOR_PROBABILITIES_KEY, MAP_ALLELE_FRACTIONS_KEY);
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        Utils.nonNull(likelihoods);

        // do not annotate the genotype fields for normal
        if (g.isHomRef()){
            return;
        }

        pi.put(NO_ARTIFACT, 0.95);
        pi.put(ART_FWD, 0.025);
        pi.put(ART_REV, 0.025);

        // We use the allele with highest LOD score
        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY, () -> null, -1);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
        final Allele altAlelle = vc.getAlternateAllele(indexOfMaxTumorLod);

        final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAlleles(g.getSampleName());
        final int numFwdAltReads = (int) bestAlleles.stream().filter(ba -> !ba.read.isReverseStrand() && ba.isInformative() && ba.allele.equals(altAlelle)).count();
        final int numRevAltReads = (int) bestAlleles.stream().filter(ba -> ba.read.isReverseStrand() && ba.isInformative() && ba.allele.equals(altAlelle)).count();
        final int numFwdReads = (int) bestAlleles.stream().filter(ba -> !ba.read.isReverseStrand() && ba.isInformative()).count();
        final int numRevReads = (int) bestAlleles.stream().filter(ba -> ba.read.isReverseStrand() && ba.isInformative()).count();
        final int numAltReads = numFwdAltReads + numRevAltReads;
        final int numReads = numFwdReads + numRevReads;

        final EnumMap<StrandArtifactZ, Double> unnormalized_posterior_probabilities = new EnumMap<>(StrandArtifactZ.class);
        final EnumMap<StrandArtifactZ, Double> maximum_a_posteriori_allele_fraction_estimates = new EnumMap<>(StrandArtifactZ.class);

        /*** Compute the posterior probability of ARTIFACT_FWD and ARTIFACT_REV; it's a double integral over f and epsilon ***/

        // the integrand is a polynomial of degree n, where n is the number of reads at the locus
        // thus to integrate exactly with Gauss-Legendre we need (n/2)+1 points
        final int numIntegPointsForAlleleFraction = numReads / 2 + 1;
        final int numIntegPointsForEpsilon = (numReads + ALPHA + BETA - 2) / 2 + 1;

        final double likelihoodForArtifactFwd = IntegrationUtils.integrate2d(
                (f,epsilon) -> getIntegrandGivenArtifact(f, epsilon, numFwdReads, numRevReads, numFwdAltReads, numRevAltReads),
                0.0, 1.0, numIntegPointsForAlleleFraction,
                0.0, 1.0, numIntegPointsForEpsilon);
        final double likelihoodForArtifactRev = IntegrationUtils.integrate2d(
                (f,epsilon) -> getIntegrandGivenArtifact(f, epsilon, numRevReads, numFwdReads, numRevAltReads, numFwdAltReads),
                0.0, 1.0, numIntegPointsForAlleleFraction,
                0.0, 1.0, numIntegPointsForEpsilon);

        unnormalized_posterior_probabilities.put(ART_FWD, pi.get(ART_FWD) * likelihoodForArtifactFwd);
        unnormalized_posterior_probabilities.put(ART_REV, pi.get(ART_REV) * likelihoodForArtifactRev);

        /*** Compute the posterior probability of NO_ARTIFACT; evaluate a single integral over the allele fraction ***/
        final double likelihoodForNoArtifact = IntegrationUtils.integrate(
                f -> getIntegrandGivenNoArtifact(f, numFwdReads, numRevReads, numFwdAltReads, numRevAltReads),
                0.0, 1.0, numIntegPointsForAlleleFraction);

        unnormalized_posterior_probabilities.put(NO_ARTIFACT, pi.get(NO_ARTIFACT) * likelihoodForNoArtifact);

        final double[] posterior_probabilities = MathUtils.normalizeFromRealSpace(
                unnormalized_posterior_probabilities.values().stream().mapToDouble(Double::doubleValue).toArray());

        /*** Compute the maximum a posteriori estimate for allele fraction given strand artifact ***/
        // For a fixed f, integrate the double integral over epsilons. This gives us the likelihood p(x^+, x^- | f, z) for a fixed f, which is proportional to
        // the posterior probability of p(f | x^+, x^-, z)
        final int numSamplePoints = 100;
        final double[] samplePoints = GATKProtectedMathUtils.createEvenlySpacedPoints(0.0, 1.0, numSamplePoints);
        double[] likelihoodsGivenForwardArtifact = new double[numSamplePoints];
        double[] likelihoodsGivenReverseArtifact = new double[numSamplePoints];

        for (int i = 0; i < samplePoints.length; i++){
            final double f = samplePoints[i];
            likelihoodsGivenForwardArtifact[i] = IntegrationUtils.integrate(
                    epsilon -> getIntegrandGivenArtifact(f, epsilon, numFwdReads, numRevReads, numFwdAltReads, numRevAltReads),
                    0.0, 1.0, numIntegPointsForEpsilon);
            likelihoodsGivenReverseArtifact[i] = IntegrationUtils.integrate(
                    epsilon -> getIntegrandGivenArtifact(f, epsilon, numRevReads, numFwdReads, numRevAltReads, numFwdAltReads),
                    0.0, 1.0, numIntegPointsForEpsilon);
        }

        final int maxAlleleFractionIndexFwd = MathUtils.maxElementIndex(likelihoodsGivenForwardArtifact);
        final int maxAlleleFractionIndexRev = MathUtils.maxElementIndex(likelihoodsGivenReverseArtifact);
        maximum_a_posteriori_allele_fraction_estimates.put(ART_FWD, samplePoints[maxAlleleFractionIndexFwd]);
        maximum_a_posteriori_allele_fraction_estimates.put(ART_REV, samplePoints[maxAlleleFractionIndexRev]);

        // In the absence of strand artifact, MAP estimate for f reduces to the sample alt allele fraction
        maximum_a_posteriori_allele_fraction_estimates.put(NO_ARTIFACT,  (double) numAltReads/numReads);

        gb.attribute(POSTERIOR_PROBABILITIES_KEY, posterior_probabilities);
        gb.attribute(MAP_ALLELE_FRACTIONS_KEY, maximum_a_posteriori_allele_fraction_estimates.values().stream().mapToDouble(Double::doubleValue).toArray());
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(POSTERIOR_PROBABILITIES_KEY, 3, VCFHeaderLineType.Float, "posterior probabilities of the presence of strand artifact"),
                new VCFFormatHeaderLine(MAP_ALLELE_FRACTIONS_KEY, 3, VCFHeaderLineType.Float, "MAP estimates of allele fraction given z"));
    }

    private double getIntegrandGivenNoArtifact(final double f, final int nPlus, final int nMinus, final int xPlus, final int xMinus){
        return MathUtils.binomialProbability(nPlus, xPlus, f) * MathUtils.binomialProbability(nMinus, xMinus, f);
    }

    private double getIntegrandGivenArtifact(final double f, final double epsilon, final int nWithArtifact, final int nNoArtifact, final int xWithArtifact, final int xNoArtifact){
        return betaPrior.density(epsilon) * MathUtils.binomialProbability(nWithArtifact, xWithArtifact, f + epsilon * (1-f)) *
                MathUtils.binomialProbability(nNoArtifact, xNoArtifact, f);
    }


}
