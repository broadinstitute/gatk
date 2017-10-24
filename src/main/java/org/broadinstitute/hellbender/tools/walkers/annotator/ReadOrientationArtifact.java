package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.readorientation.Hyperparameters;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.State;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;

/**
 * Created by tsato on 10/20/17.
 */
public class ReadOrientationArtifact extends GenotypeAnnotation implements StandardMutectAnnotation {
    private File hyperparameterTable;

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.READ_ORIENTATION_POSTERIOR_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(GATKVCFConstants.READ_ORIENTATION_POSTERIOR_KEY, 2, 
                VCFHeaderLineType.Float, "posterior probabilities of read orientation artifact [p(f1r2), p(f1r1)]"));
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

        // do not annotate the genotype fields for the normal sample
        if (g.isHomRef()){
            return;
        }

        // skip everything but a SNP
        if (! vc.isSNP()){
            return;
        }

        // must get a depth - how shall I do it?

        final int altAlleleIndex = 1;
        // TODO: extract this into a method - getKmerAround() is lacking
        ref.setWindow(1, 1);
        final String refContext = new String(ref.getBases()); // this might take longer we think
        Utils.validate(refContext.length() == 3,"kmer must have length 3");

        if (refContext.contains("N")) {
            return;
        }

        final Nucleotide refBase = Nucleotide.valueOf(refContext.substring(1,2));

        // Using a stream over a collection of a parameterized type is messy - use a for loop instead
        final int[] baseCounts = new int[4];
        final int[] f1r2Counts = new int[4];
        final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAlleles(g.getSampleName());
        for (ReadLikelihoods<Allele>.BestAllele bestAllele : bestAlleles){
            final Allele allele = bestAllele.allele;
            Nucleotide base = allele.isReference() ? refBase : Nucleotide.valueOf(allele.getBaseString());;
            baseCounts[base.ordinal()] += 1;
            if (! ReadUtils.isF2R1(bestAllele.read)){
                f1r2Counts[base.ordinal()] += 1;
            }
        }

        final List<Hyperparameters> hyperparametersForReadOrientaitonModel =
                Hyperparameters.readHyperparameters(hyperparameterTable);

        final Nucleotide refAllele = refBase;
        final Nucleotide altAllele = Nucleotide.valueOf(vc.getAlternateAllele(0).toString());

        final State[] artifactStatesOfInterest = State.getArtifactStatesOfInterest(altAllele);

        final Optional<Hyperparameters> optionalHyps = hyperparametersForReadOrientaitonModel.stream()
                .filter(h -> h.getReferenceContext().equals(refContext))
                .findFirst();

        if (! optionalHyps.isPresent()){
            // without the hyperparameters for this particular reference context we cannot compute the posterior probabilities
            return;
        }

        Hyperparameters hyps = optionalHyps.get();

        // A by K matrix of prior probabilities over K latent states, given allele a \in A
        // \pi_{ak} is the prior probability of state k given observed allele a.
        final double[] pi = hyps.getPi();

        // a vector of length K, the probability of drawing an alt read (i.e. allele fraction) given z
        final double[] f = hyps.getF();

        // a vector of length K, the probability of drawing an F1R2 alt read given z
        final double[] theta = hyps.getTheta();

        final double[] log10UnnormalizedPosteriorProbabilities = new double[LearnHyperparametersEngine.NUM_STATES];

        // gather data
        final int depth = (int) MathUtils.sum(g.getAD());
        final int altDepth = g.getAD()[altAlleleIndex];
        final int f1r2AltCount = f1r2Counts[altAllele.ordinal()];

        List<State> impossibleStates = State.getImpossibleStates(refAllele);
        for (State z : State.values()){
            final int k = z.ordinal();
            if (impossibleStates.contains(z)){
                log10UnnormalizedPosteriorProbabilities[k] = Double.NEGATIVE_INFINITY;
                continue;
            }

            final int m_nk; // alt depth for state k
            final int x_nk; // alt f1r2 depth for state k
            final int r = depth; // depth at this site

            if (State.getNonArtifactStates().contains(z) || State.getAltAlleleOfTransition(z) == altAllele){
                // We are in a non-artifact state {germline het, somatic het, hom ref, hom var}
                // Or we are in an artifact state that matches the alternate allele
                // e.g. z = F1R2_c when alt allele = C
                m_nk = altDepth;
                x_nk = f1r2AltCount;
            } else {
                // We approximate the counts of bases that are neither ref nor alt (if ref = A, alt = C, then {G, T})
                // Although it would be better to use the exact count of each base, this would require
                // making a new annotation and cluttering the vcf, and I doube that that will improve the results
                m_nk = baseCounts[State.getAltAlleleOfTransition(z).ordinal()];
                x_nk = f1r2Counts[State.getAltAlleleOfTransition(z).ordinal()];
            }

            log10UnnormalizedPosteriorProbabilities[k] =
                    LearnHyperparametersEngine.computePosterior(z, m_nk, x_nk, r, pi[k], f[k], theta[k]);
        }

        final double[] posteriorProbabilties = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedPosteriorProbabilities);

        final double f1r2PosteriorProbability = posteriorProbabilties[artifactStatesOfInterest[0].ordinal()];
        final double f2r1PosteriorProbability = posteriorProbabilties[artifactStatesOfInterest[1].ordinal()];

        gb.attribute(GATKVCFConstants.READ_ORIENTATION_POSTERIOR_KEY,
                new double[]{f1r2PosteriorProbability, f2r1PosteriorProbability});
    }

    /** Part of the hack scheme that is required until we can pass an argument directly ot an annotation class **/
    public void setHyperparameterTable(final File hyperparameterTable){
        this.hyperparameterTable = hyperparameterTable;
    }

}
