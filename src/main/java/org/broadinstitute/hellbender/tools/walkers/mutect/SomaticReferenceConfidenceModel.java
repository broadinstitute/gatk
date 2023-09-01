package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceResult;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class SomaticReferenceConfidenceModel extends ReferenceConfidenceModel {

    private final SampleList samples;
    private final Optional<BetaDistributionShape> afPrior;



    /**
     * Create a new ReferenceConfidenceModel
     *  @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     * @param minAF soft threshold for allele fractions -- above this value prior is nearly flat, below, prior is nearly zero
     * @param refModelDelQual reference model deletion quality (if constant)
     * @param useSoftClippedBases should soft clipped bases be counted against the reference
     * @param isFlowBasedModel is the error model flow based
     */
    SomaticReferenceConfidenceModel(final SampleList samples, final SAMFileHeader header, final int indelInformativeDepthIndelSize,
                                    final double minAF, final byte refModelDelQual, final boolean useSoftClippedBases,
                                    final boolean isFlowBasedModel){
        super(samples, header, indelInformativeDepthIndelSize, 0, refModelDelQual, useSoftClippedBases, isFlowBasedModel);
        Utils.validateArg(minAF >= 0.0 && minAF < 1, "minAF must be < 1 and >= 0");

        // To softly cut off allele fractions below minAF, we use a Beta prior of the form Beta(1+epsilon, 1); that is
        // the prior on allele fraction f is proportional to f^epsilon.  If epsilon is small this prior vanishes as f -> 0
        // and very rapidly becomes flat.  We choose epsilon such that minAF^epsilon = 0.5.
        afPrior = minAF == 0.0 ? Optional.empty() : Optional.of(new BetaDistributionShape(1 - Math.log(2)/Math.log(minAF), 1));
        this.samples = samples;
    }

        /**
         * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
         *
         * @param ploidy target sample ploidy.
         * @param pileup the read backed pileup containing the data we want to evaluate
         * @param refBase the reference base at this pileup position
         * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
         * @param hqSoftClips running average data structure (can be null) to collect information about the number of high quality soft clips
         * @return a RefVsAnyResult genotype call.
         */
    @Override
    public ReferenceConfidenceResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
                                                                       final ReadPileup pileup,
                                                                       final byte refBase,
                                                                       final byte minBaseQual,
                                                                       final MathUtils.RunningAverage hqSoftClips,
                                                                       final boolean readsWereRealigned) {

        final SomaticRefVsAnyResult result = new SomaticRefVsAnyResult();
        final Map<String, List<GATKRead>> perSampleReadMap = new HashMap<>();
        perSampleReadMap.put(samples.getSample(0), pileup.getReads());

        final List<Byte> altQuals = new ArrayList<>(pileup.size() / 20);

        for (final PileupElement element : pileup) {
            if (!element.isDeletion() && element.getQual() <= minBaseQual) {
                continue;
            }

            final boolean isAlt = readsWereRealigned ? isAltAfterAssembly(element, refBase) : isAltBeforeAssembly(element, refBase);
            if (isAlt) {
                altQuals.add(element.getQual());
                result.nonRefDepth++;
            } else {
                result.refDepth++;
            }
        }

        final double logOdds = Mutect2Engine.logLikelihoodRatio(result.refDepth, altQuals, 1, afPrior);
        result.lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        result.lods.set(Allele.NON_REF_ALLELE, logOdds);

        return result;
    }

    @Override
    public void addGenotypeData(final ReferenceConfidenceResult result, final GenotypeBuilder gb) {
        gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, MathUtils.logToLog10(((SomaticRefVsAnyResult)result).lods.get(Allele.NON_REF_ALLELE)));
    }

    @Override
    public void doIndelRefConfCalc(final int ploidy, final byte[] ref, final ReadPileup pileup, final int refOffset, final ReferenceConfidenceResult homRefCalc) {
        //NOTE:
        // For germline we evaluate the alternative indel reference confidence model, compare with the SNP ref conf
        // model results, and return the less confident likelihoods. The existing indel model finds the number of reads
        // spanning the current position that are informative for indels of size [-10, 10]bp and calculates a diploid
        // genotype likelihood with a constant indel "quality" and up to 40 informative reads. I'm not convinced that
        // low allele fraction variants or errors derived from PCR error are well represented in that model.
        // Additionally, the first application of this somatic ref conf is for mitochondrial calling. The MT reference
        // is quite high complexity so the SNP model should dominate. Until we develop a better model for somatic
        // indels, we will rely on the SNP model for all applications and this method (which is called by the parent class)
        // will be a noop.
    }

}
