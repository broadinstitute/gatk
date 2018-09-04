package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.annotator.StrandArtifact.ArtifactState.*;

/**
 * Annotations for strand artifact filter calculated by a Bayesian model described in <a href='https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf' target='_blank'>https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf</a>.
 *
 * <p>Output consists of two three-element arrays.  The elements of SA_POST_PROB are, respectively, the posterior probabilities that there exists an artifact on the forward strand, an artifact on the reverse strand, or no artifact.
 * These probabilities are normalized to sum to 1.  The elements of SA_MAP_AF are the max a posteriori estimates of the variant allele fraction given a
 * forward strand artifact, a reverse strand artifact, and no artifact, respectively.  For example SA_POST_PROB = 0.9, 0.02, 0.08 and
 * SA_MAP_AF = 0.01, 0.1, 0.2 together mean that the apparent variant is most likely a forward strand artifact, but if it is actually a true
 * variant its allele fraction is most likely 0.2.</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Annotations for strand artifact filter (SA_POST_PROB, SA_MAP_AF)")
public class StrandArtifact extends GenotypeAnnotation implements StandardMutectAnnotation {
    protected final OneShotLogger warning = new OneShotLogger(this.getClass());

    public static final String POSTERIOR_PROBABILITIES_KEY = "SA_POST_PROB";
    public static final String MAP_ALLELE_FRACTIONS_KEY = "SA_MAP_AF";

    // pseudocounts for the beta distribution over epsilon
    // alpha > 0 and beta > 0. alpha = beta = 1 gives us the flat prior
    // give more prior weight to beta (i.e. pseudocount for tails) so that the peak shifts towards 0
    private static final int ALPHA = 1;
    private static final int BETA = 6;
    private static final BetaDistribution betaPrior = new BetaDistribution(null, ALPHA, BETA);

    private static final double PRIOR_PROBABILITY_OF_NO_ARTIFACT = 0.95;
    private static final double PRIOR_PROBABILITY_OF_ARTIFACT = (1 - PRIOR_PROBABILITY_OF_NO_ARTIFACT)/2;

    // Apache Commons uses naive recursion for Gauss-Legendre and is prone to stack overflows
    // capping the number of subdivisions is a stop-gap for a more principled integration scheme
    private static final int MAX_GAUSS_LEGENDRE_POINTS = 100;

    // the discrete latent variable that represents the artifact state
    public enum ArtifactState {
        ART_FWD, ART_REV, NO_ARTIFACT
    }

    private static final EnumMap<ArtifactState, Double> pi = new EnumMap<>(ArtifactState.class);

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

        // do not annotate the genotype of the normal sample
        if (g.isHomRef()){
            return;
        }

        pi.put(NO_ARTIFACT, PRIOR_PROBABILITY_OF_NO_ARTIFACT);
        pi.put(ART_FWD, PRIOR_PROBABILITY_OF_ARTIFACT);
        pi.put(ART_REV, PRIOR_PROBABILITY_OF_ARTIFACT);

        // We use the allele with highest LOD score
        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY, () -> null, -1);

        if (tumorLods==null) {
            warning.warn("One or more variant contexts is missing the 'TLOD' annotation, StrandArtifact will not be computed for these VariantContexts");
            return;
        }
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
        final Allele altAlelle = vc.getAlternateAllele(indexOfMaxTumorLod);

        final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(g.getSampleName());
        final int numFwdAltReads = (int) bestAlleles.stream().filter(ba -> !ba.read.isReverseStrand() && ba.isInformative() && ba.allele.equals(altAlelle)).count();
        final int numRevAltReads = (int) bestAlleles.stream().filter(ba -> ba.read.isReverseStrand() && ba.isInformative() && ba.allele.equals(altAlelle)).count();
        final int numFwdReads = (int) bestAlleles.stream().filter(ba -> !ba.read.isReverseStrand() && ba.isInformative()).count();
        final int numRevReads = (int) bestAlleles.stream().filter(ba -> ba.read.isReverseStrand() && ba.isInformative()).count();
        final int numAltReads = numFwdAltReads + numRevAltReads;
        final int numReads = numFwdReads + numRevReads;

        final EnumMap<ArtifactState, Double> unnormalizedPosterior = new EnumMap<>(ArtifactState.class);
        final EnumMap<ArtifactState, Double> estimatedAlleleFractions = new EnumMap<>(ArtifactState.class);

        /** Compute the posterior probabilities of ARTIFACT_FWD and ARTIFACT_REV states, which entails
         *  a double integral over the latent variables f and epsilon
         **/

        // The integrand is a polynomial of degree n, where n is the number of reads at the locus.
        // To integrate exactly with Gauss-Legendre, we need (n/2)+1 points
        final int numIntegPointsForAlleleFraction = Math.min(numReads / 2 + 1, MAX_GAUSS_LEGENDRE_POINTS);
        final int numIntegPointsForEpsilon = Math.min((numReads + ALPHA + BETA - 2) / 2 + 1, MAX_GAUSS_LEGENDRE_POINTS);

        final double likelihoodForArtifactFwd = IntegrationUtils.integrate2d(
                (f,epsilon) -> getIntegrandGivenArtifact(f, epsilon, numFwdReads, numRevReads, numFwdAltReads, numRevAltReads),
                0.0, 1.0, numIntegPointsForAlleleFraction,
                0.0, 1.0, numIntegPointsForEpsilon);
        final double likelihoodForArtifactRev = IntegrationUtils.integrate2d(
                (f,epsilon) -> getIntegrandGivenArtifact(f, epsilon, numRevReads, numFwdReads, numRevAltReads, numFwdAltReads),
                0.0, 1.0, numIntegPointsForAlleleFraction,
                0.0, 1.0, numIntegPointsForEpsilon);

        unnormalizedPosterior.put(ART_FWD, pi.get(ART_FWD) * likelihoodForArtifactFwd);
        unnormalizedPosterior.put(ART_REV, pi.get(ART_REV) * likelihoodForArtifactRev);

        /*** Compute the posterior probability of NO_ARTIFACT; evaluate a single integral over the allele fraction ***/
        final double likelihoodForNoArtifact = IntegrationUtils.integrate(
                f -> getIntegrandGivenNoArtifact(f, numFwdReads, numRevReads, numFwdAltReads, numRevAltReads),
                0.0, 1.0, numIntegPointsForAlleleFraction);

        unnormalizedPosterior.put(NO_ARTIFACT, pi.get(NO_ARTIFACT) * likelihoodForNoArtifact);

        final double[] posteriorProbabilities = MathUtils.normalizeFromRealSpace(
                unnormalizedPosterior.values().stream().mapToDouble(Double::doubleValue).toArray());

        /*** Compute the maximum a posteriori estimate for allele fraction given strand artifact ***/
        // For a fixed f, integrate the double integral over epsilons. This gives us the likelihood p(x^+, x^- | f, z) for a fixed f, which is proportional to
        // the posterior probability of p(f | x^+, x^-, z)
        final int numSamplePoints = 100;
        final double[] samplePoints = MathUtils.createEvenlySpacedPoints(0.0, 1.0, numSamplePoints);
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
        estimatedAlleleFractions.put(ART_FWD, samplePoints[maxAlleleFractionIndexFwd]);
        estimatedAlleleFractions.put(ART_REV, samplePoints[maxAlleleFractionIndexRev]);

        // In the absence of strand artifact, MAP estimate for f reduces to the sample alt allele fraction
        estimatedAlleleFractions.put(NO_ARTIFACT,  (double) numAltReads/numReads);

        gb.attribute(POSTERIOR_PROBABILITIES_KEY, posteriorProbabilities);
        gb.attribute(MAP_ALLELE_FRACTIONS_KEY, estimatedAlleleFractions.values().stream().mapToDouble(Double::doubleValue).toArray());
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
