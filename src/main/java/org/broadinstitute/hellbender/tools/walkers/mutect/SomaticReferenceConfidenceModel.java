package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RefVsAnyResult;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceResult;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class SomaticReferenceConfidenceModel extends ReferenceConfidenceModel {

    private final SampleList samples;
    private final SomaticGenotypingEngine genotypingEngine;



    /**
     * Create a new ReferenceConfidenceModel
     *
     * @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     */
    SomaticReferenceConfidenceModel(final SampleList samples,
                                    final SAMFileHeader header,
                                    final int indelInformativeDepthIndelSize,
                                    final SomaticGenotypingEngine genotypingEngine){
        super(samples, header, indelInformativeDepthIndelSize, 0);
        this.samples = samples;
        this.genotypingEngine = genotypingEngine;
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
        final ReadLikelihoods<Allele> readLikelihoods = new ReadLikelihoods<>(samples, new IndexedAlleleList<>(Arrays.asList(Allele.create(refBase,true), Allele.NON_REF_ALLELE)), perSampleReadMap);
        final ReadLikelihoods<Allele> readLikelihoods2 = new ReadLikelihoods<>(samples, new IndexedAlleleList<>(Arrays.asList(Allele.create(refBase,true), Allele.NON_REF_ALLELE)), perSampleReadMap);
        final Iterator<PileupElement> pileupIter = pileup.iterator();
        for (int i = 0; i < pileup.size(); i++) {
            final PileupElement element = pileupIter.next();
            if (!element.isDeletion() && element.getQual() <= minBaseQual) {
                continue;
            }
            final boolean isAlt = readsWereRealigned ? isAltAfterAssembly(element, refBase) : isAltBeforeAssembly(element, refBase);
            final double nonRefLikelihood;
            final double refLikelihood;
            if (isAlt) {
                nonRefLikelihood = QualityUtils.qualToProbLog10(element.getQual());
                refLikelihood = QualityUtils.qualToErrorProbLog10(element.getQual()) + MathUtils.LOG10_ONE_THIRD;
                result.nonRefDepth++;
            } else {
                nonRefLikelihood = QualityUtils.qualToErrorProbLog10(element.getQual()) + MathUtils.LOG10_ONE_THIRD;
                refLikelihood = QualityUtils.qualToProbLog10(element.getQual());
                result.refDepth++;
            }
            readLikelihoods.sampleMatrix(0).set(0, i, nonRefLikelihood);
            readLikelihoods2.sampleMatrix(0).set(0, i, refLikelihood);
            readLikelihoods2.sampleMatrix(0).set(1, i, nonRefLikelihood);
        }
        result.lods = genotypingEngine.somaticLog10Odds(readLikelihoods.sampleMatrix(0));
        PerAlleleCollection<Double> lods2 = genotypingEngine.somaticLog10Odds(readLikelihoods2.sampleMatrix(0));
        result.lods = lods2;
        return result;
    }

    @Override
    public void addGenotypeData(final ReferenceConfidenceResult result, final GenotypeBuilder gb) {
        gb.attribute(GATKVCFConstants.TUMOR_LOD_KEY, ((SomaticRefVsAnyResult)result).lods.get(Allele.NON_REF_ALLELE));
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
