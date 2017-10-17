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
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.ArtifactState;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.CollectDataForReadOrientationFilter.REF_CONTEXT_PADDING_ON_EACH_SIDE;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.CollectDataForReadOrientationFilter.REGULAR_BASES;

/**
 * Created by tsato on 10/20/17.
 */
public class ReadOrientationArtifact extends GenotypeAnnotation implements StandardMutectAnnotation {
    private File hyperparameterTable;
    public static int INDEX_OF_F1R2_ARTIFACT = 0;
    public static int INDEX_OF_F2R1_ARTIFACT = 1;


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

        if (hyperparameterTable == null){
            return;
        }

        // must get a depth - how shall I do it?

        final int altAlleleIndex = 1;
        final String refContext = ReferenceContext.extractKmer(vc.getStart(), ref, REF_CONTEXT_PADDING_ON_EACH_SIDE);
        Utils.validate(refContext.length() == 2*REF_CONTEXT_PADDING_ON_EACH_SIDE + 1,
                String.format("kmer must have length %d but got %d", 2* REF_CONTEXT_PADDING_ON_EACH_SIDE + 1, refContext.length()));

        if (refContext.contains("N")) {
            return;
        }

        final Nucleotide refAllele = Nucleotide.valueOf(refContext.substring(REF_CONTEXT_PADDING_ON_EACH_SIDE,REF_CONTEXT_PADDING_ON_EACH_SIDE+1));
        Utils.validate(refAllele == Nucleotide.valueOf(vc.getReference().toString().replace("*", "")),
                String.format("ref allele in the kmer, %s, does not match the ref allele in the variant context, %s",
                        refAllele, vc.getReference().toString().replace("*", "")));

        // Using a stream over a collection of a parameterized type is messy - use a for loop instead
        final int[] baseCounts = new int[REGULAR_BASES.size()];
        final int[] f1r2Counts = new int[REGULAR_BASES.size()];
        final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAlleles(g.getSampleName());
        for (ReadLikelihoods<Allele>.BestAllele bestAllele : bestAlleles){
            final Allele allele = bestAllele.allele;

            if (!bestAllele.isInformative() || allele.length() != 1){
                continue;
            }

            final GATKRead read = bestAllele.read;

            // We will throw out any read that does not overlap the locus (this happens, presumably, when there are other reads
            // that are in phase with this particular read and the engine imputes the missing base according to their haplotype)
            if (read.getEnd() < vc.getStart()){
                continue;
            }

            // final int offset = vc.getStart() - read.getStart();
            final AlignmentStateMachine asm = new AlignmentStateMachine(read);
            while ( asm.stepForwardOnGenome() != null && asm.getGenomePosition() < vc.getStart()) {

            }

            final int minimumBaseQuality = 20;
            final int readOffset = asm.getReadOffset();
            if (asm.getGenomePosition() != vc.getStart()){
                System.out.println(String.format("asm.getGenomePosition() %s:%d != vc.getStart() %d\nread = %s",
                        asm.getContig(), asm.getGenomePosition(), vc.getStart(),read.getName()));
            }

            if (asm.getGenomePosition() == vc.getStart() && read.getBaseQuality(readOffset) < minimumBaseQuality){
                continue;
            }

            Nucleotide base = allele.isReference() ? refAllele : Nucleotide.valueOf(allele.getBaseString());;
            baseCounts[base.ordinal()] += 1;
            if (! ReadUtils.isF2R1(bestAllele.read)){
                f1r2Counts[base.ordinal()] += 1;
            }
        }

        final List<Hyperparameters> hyperparametersForReadOrientaitonModel =
                Hyperparameters.readHyperparameters(hyperparameterTable);

        final Nucleotide altAllele = Nucleotide.valueOf(vc.getAlternateAllele(0).toString());

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

        // gather data
        final int depth = (int) MathUtils.sum(baseCounts);
        final int altDepth = baseCounts[altAllele.ordinal()];
        final int f1r2AltCount = f1r2Counts[altAllele.ordinal()];

        final double[] log10UnnormalizedPosteriorProbabilities = LearnHyperparametersEngine.computeLog10Responsibilities(
                refAllele, altAllele, altDepth, f1r2AltCount, depth, pi);

        // we want the posterior of artifacts given that the site is not hom ref
        log10UnnormalizedPosteriorProbabilities[LearnHyperparametersEngine.ArtifactState.HOM_REF.ordinal()] = Double.NEGATIVE_INFINITY;

        final double[] posteriorProbabilities = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedPosteriorProbabilities);

        final double[] posteriorProbabilitiesOfArtifact = new double[2];
        posteriorProbabilitiesOfArtifact[INDEX_OF_F1R2_ARTIFACT] = posteriorProbabilities[ArtifactState.getF1R2StateOfInterest(altAllele).ordinal()];
        posteriorProbabilitiesOfArtifact[INDEX_OF_F2R1_ARTIFACT] = posteriorProbabilities[ArtifactState.getF2R1StateOfInterest(altAllele).ordinal()];

        // TODO: change to one key saying F1R2, and the other probability
        gb.attribute(GATKVCFConstants.READ_ORIENTATION_POSTERIOR_KEY, posteriorProbabilitiesOfArtifact);
    }

    /** Part of the hack scheme that is required until we can pass an argument directly ot an annotation class **/
    public void setHyperparameterTable(final File hyperparameterTable){
        this.hyperparameterTable = hyperparameterTable;
    }

}
