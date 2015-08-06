package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.indels.PairHMMIndelErrorModel;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.genotyper.ProbabilityVector;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 * This is a site based implementation of an Error Model. The error model is a probability
 * distribution for the site given the phred scaled quality.
 */
public class ErrorModel {
    private byte maxQualityScore;
    private byte minQualityScore;
    private byte phredScaledPrior;
    private double log10minPower;
    private int refDepth;
    private boolean hasData = false;
    private ProbabilityVector probabilityVector;
    private static final boolean compressRange = false;
    
    private static final double log10MinusE = Math.log10(Math.exp(1.0));
    private static final boolean DEBUG = false;
    /**
     * Calculates the probability of the data (reference sample reads) given the phred scaled site quality score.
     * 
     * @param UAC                           Argument Collection
     * @param refSamplePileup            Reference sample pileup
     * @param refSampleVC                VC with True alleles in reference sample pileup
     */
    public ErrorModel(final UnifiedArgumentCollection UAC,
                      final ReadPileup refSamplePileup,
                      final VariantContext refSampleVC,
                      final ReferenceContext refContext) {
        this.maxQualityScore = UAC.maxQualityScore;
        this.minQualityScore = UAC.minQualityScore;
        this.phredScaledPrior = UAC.phredScaledPrior;
        log10minPower = Math.log10(UAC.minPower);

        PairHMMIndelErrorModel pairModel = null;
        LinkedHashMap<Allele, Haplotype> haplotypeMap = null;
        double[][] perReadLikelihoods = null;

        final double[] model = new double[maxQualityScore+1];
        Arrays.fill(model, Double.NEGATIVE_INFINITY);

        boolean hasCalledAlleles = false;

        final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        if (refSampleVC != null) {

            for (final Allele allele : refSampleVC.getAlleles()) {
                if (allele.isCalled()) {
                    hasCalledAlleles = true;
                    break;
                }
            }
            haplotypeMap = new LinkedHashMap<>();
            if (refSampleVC.isIndel()) {
                pairModel = new PairHMMIndelErrorModel(UAC.INDEL_GAP_OPEN_PENALTY, UAC.INDEL_GAP_CONTINUATION_PENALTY,
                        UAC.OUTPUT_DEBUG_INDEL_INFO, UAC.pairHMM);
                IndelGenotypeLikelihoodsCalculationModel.getHaplotypeMapFromAlleles(refSampleVC.getAlleles(), refContext, refContext.getInterval(), haplotypeMap); // will update haplotypeMap adding elements
            }
        }

        final double p = QualityUtils.qualToErrorProbLog10((byte) (maxQualityScore - minQualityScore));
        if (refSamplePileup == null || refSampleVC == null  || !hasCalledAlleles) {
            for (byte q=minQualityScore; q<=maxQualityScore; q++) {
                // maximum uncertainty if there's no ref data at site
                model[q] = p;
            }
            this.refDepth = 0;
        }
        else {
            hasData = true;
            int matches = 0;
            int coverage = 0;

            final Allele refAllele = refSampleVC.getReference();

            if ( refSampleVC.isIndel()) {
                //perReadLikelihoods = new double[readCounts.length][refSampleVC.getAlleles().size()];
                final int eventLength = IndelGenotypeLikelihoodsCalculationModel.getEventLength(refSampleVC.getAlleles());
                if (!haplotypeMap.isEmpty()) {
                    perReadLikelihoods = pairModel.computeGeneralReadHaplotypeLikelihoods(refSamplePileup, haplotypeMap, refContext, eventLength, perReadAlleleLikelihoodMap);
                }
            }
            int idx = 0;
            for (final PileupElement refPileupElement : refSamplePileup) {
                if (DEBUG) {
                    System.out.println(refPileupElement.toString());
                }
                boolean isMatch = false;
                for (final Allele allele : refSampleVC.getAlleles()) {
                    final boolean m = pileupElementMatches(refPileupElement, allele, refAllele, refContext.getBase());
                    if (DEBUG) {
                        System.out.println(m);
                    }
                    isMatch |= m;
                }
                if (refSampleVC.isIndel() && !haplotypeMap.isEmpty()) {
                    // ignore match/mismatch if reads, as determined by their likelihood, are not informative
                    final double[] perAlleleLikelihoods = perReadLikelihoods[idx++];
                    if (!isInformativeElement(perAlleleLikelihoods)) {
                        matches++;
                    } else {
                        matches += (isMatch ? 1 : 0);
                    }

                }   else {
                    matches += (isMatch?1:0);
                }
                coverage++;
            }

            final int mismatches = coverage - matches;
            //System.out.format("Cov:%d match:%d mismatch:%d\n",coverage, matches, mismatches);
            for (byte q=minQualityScore; q<=maxQualityScore; q++) {
                if (coverage==0) {
                    model[q] = p;
                } else {
                    model[q] = log10PoissonProbabilitySiteGivenQual(q, coverage, mismatches);
                }
            }
            this.refDepth = coverage;
        }
        
        // compress probability vector
        this.probabilityVector = new ProbabilityVector(model, compressRange);
    }

    private boolean isInformativeElement(final double[] likelihoods) {
        // if likelihoods are the same, they're not informative
        final double thresh = 0.1;
        final int maxIdx = MathUtils.maxElementIndex(likelihoods);
        final int minIdx = MathUtils.minElementIndex(likelihoods);
        if (likelihoods[maxIdx]-likelihoods[minIdx]< thresh) {
            return false;
        } else {
            return true;
        }
    }
    /**
     * Simple constructor that just takes a given log-probability vector as error model.
     * Only intended for unit testing, not general usage.
     * @param pvector       Given vector of log-probabilities
     *
     */
    public ErrorModel(final double[] pvector) {
        this.maxQualityScore = (byte)(pvector.length-1);
        this.minQualityScore = 0;
        this.probabilityVector = new ProbabilityVector(pvector, compressRange);
        this.hasData = true;

    }

    public static boolean pileupElementMatches(final PileupElement pileupElement, final Allele allele, final Allele refAllele, final byte refBase) {
        if (DEBUG) {
            System.out.format("PE: base:%s isNextToDel:%b isNextToIns:%b eventBases:%s eventLength:%d Allele:%s RefAllele:%s\n",
                    pileupElement.getBase(), pileupElement.isBeforeDeletionStart(),
                    pileupElement.isBeforeInsertion(), pileupElement.getBasesOfImmediatelyFollowingInsertion(), pileupElement.getLengthOfImmediatelyFollowingIndel(), allele.toString(), refAllele.toString());
        }

        //pileupElement.
        // if test allele is ref, any base mismatch, or any insertion/deletion at start of pileup count as mismatch
        if (allele.isReference()) {
            // for a ref allele, any base mismatch or new indel is a mismatch.
            if(allele.getBases().length>0)
                // todo - can't check vs. allele because allele is not padded so it doesn't include the reference base at this location
                // could clean up/simplify this when unpadding is removed
            {
                return (pileupElement.getBase() == refBase && !pileupElement.isBeforeInsertion() && !pileupElement.isBeforeDeletionStart());
            } else
                // either null allele to compare, or ref/alt lengths are different (indel by definition).
                // if we have an indel that we are comparing against a REF allele, any indel presence (of any length/content) is a mismatch
            {
                return (!pileupElement.isBeforeInsertion() && !pileupElement.isBeforeDeletionStart());
            }
        }

        // for non-ref alleles to compare:
        if (refAllele.getBases().length == allele.getBases().length)
            // alleles have the same length (eg snp or mnp)
        {
            return pileupElement.getBase() == allele.getBases()[0];
        }

        // for non-ref alleles,
        final byte[] alleleBases = allele.getBases();
        final int eventLength = alleleBases.length - refAllele.getBases().length;
        if (eventLength < 0 && pileupElement.isBeforeDeletionStart() && pileupElement.getLengthOfImmediatelyFollowingIndel() == -eventLength) {
            return true;
        }

                if (eventLength > 0 && pileupElement.isBeforeInsertion() &&
                Arrays.equals(pileupElement.getBasesOfImmediatelyFollowingInsertion().getBytes(), Arrays.copyOfRange(alleleBases, 1, alleleBases.length))) // allele contains ref byte, but pileupElement's event bases doesn't
                {
                    return true;
                }

        return false;
    }


    /**
     * What's the log-likelihood that a site's quality is equal to q? If we see N observations and n mismatches,
     * and assuming each match is independent of each other and that the match probability is just dependent of
     * the site quality, so p = 10.^-q/10.
     * Since we'll normally have relatively high Q sites and deep coverage in reference samples (ie p small, N high),
     * to avoid underflows we'll use the Poisson approximation with lambda = N*p.
     * Hence, the log-likelihood of q i.e. Pr(Nmismatches = n | SiteQ = q) ~ Poisson(n | lambda = p*N) with p as above.
     * @param q                     Desired q to get likelihood from
     * @param coverage              Total coverage
     * @param mismatches            Number of mismatches
     * @return                      Likelihood of observations as a function of q
     */
    private double log10PoissonProbabilitySiteGivenQual(final byte q, final int coverage, final int mismatches) {
        // same as   log10ProbabilitySiteGivenQual but with Poisson approximation to avoid numerical underflows
        final double lambda = QualityUtils.qualToErrorProb(q) * (double )coverage;
        // log10(e^-lambda*lambda^k/k!) = -lambda + k*log10(lambda) - log10factorial(k)
        return Math.log10(lambda)*mismatches - lambda*log10MinusE- MathUtils.log10Factorial(mismatches);
    }

    public double getSiteLogErrorProbabilityGivenQual (final int qual) {
        return probabilityVector.getLogProbabilityForIndex(qual);
    }

    public byte getMaxQualityScore() {
        return maxQualityScore;
    }

    public byte getMinQualityScore() {
        return minQualityScore;
    }

    public int getMinSignificantQualityScore() {
        return new ProbabilityVector(probabilityVector,true).getMinVal();
    }

    public int getMaxSignificantQualityScore() {
        return new ProbabilityVector(probabilityVector,true).getMaxVal();
    }

    public int getReferenceDepth() {
        return refDepth;
    }
    public boolean hasData() {
        return hasData;
    }

    public ProbabilityVector getErrorModelVector() {
        return probabilityVector;
    }

    public String toString() {
        final StringBuilder result = new StringBuilder("(");
        boolean skipComma = true;
        for (final double v : probabilityVector.getProbabilityVector()) {
            if (skipComma) {
                skipComma = false;
            }
            else {
                result.append(",");
            }
            result.append(String.format("%.4f", v));
        }
        result.append(")");
        return result.toString();
    }
    
    public static int getTotalReferenceDepth(final HashMap<String, ErrorModel> perLaneErrorModels) {
        int n=0;
        for (final ErrorModel e : perLaneErrorModels.values()) {
            n += e.getReferenceDepth();
        }
        return n;
    }

}
