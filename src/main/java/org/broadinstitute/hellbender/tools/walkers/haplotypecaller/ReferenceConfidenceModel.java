package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel;
import org.broadinstitute.hellbender.tools.walkers.variantutils.PosteriorProbabilitiesUtils;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.*;

/**
 * Code for estimating the reference confidence
 *
 * This code can estimate the probability that the data for a single sample is consistent with a
 * well-determined REF/REF diploid genotype.
 *
 */
public class ReferenceConfidenceModel {

    // Annotation used to cache reference confidence information
    public static final String INDEL_INFORMATIVE_BASES_CACHE_ATTRIBUTE_NAME = "IDL";
    public static final boolean USE_CACHED_READ_INDEL_INFORMATIVENESS_VALUES = true;
    private final SampleList samples;
    private final int indelInformativeDepthIndelSize;
    private final int numRefSamplesForPrior;

    private final PosteriorProbabilitiesUtils.PosteriorProbabilitiesOptions options;

    /**
     * Surrogate quality score for no base calls.
     * <p>
     * This is the quality assigned to deletion (so without its own base-call quality) pile-up elements,
     * when assessing the confidence on the hom-ref call at that site.
     * </p>
     */
    private static final byte REF_MODEL_DELETION_QUAL = 30;

    /**
     * Base calls with quality threshold lower than this number won't be considered when assessing the
     * confidence on the hom-ref call.
     */
    private static final byte BASE_QUAL_THRESHOLD = 6;

    /**
     * Only base calls with quality strictly greater than this constant,
     * will be considered high quality if they are part of a soft-clip.
     */
    private static final byte HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD = 28;

    //TODO change this: https://github.com/broadinstitute/gsa-unstable/issues/1108
    protected static final int MAX_N_INDEL_INFORMATIVE_READS = 40; // more than this is overkill because GQs are capped at 99 anyway

    private static final int INITIAL_INDEL_LK_CACHE_PLOIDY_CAPACITY = 20;
    private static GenotypeLikelihoods[][] indelPLCache = new GenotypeLikelihoods[INITIAL_INDEL_LK_CACHE_PLOIDY_CAPACITY + 1][];

    /**
     * Indel error rate for the indel model used to assess the confidence on the hom-ref call.
     */
    private static final double INDEL_ERROR_RATE = -4.5; // 10^-4.5 indel errors per bp

    /**
     * Phred scaled qual value that corresponds to the {@link #INDEL_ERROR_RATE indel error rate}.
     */
    private static final byte INDEL_QUAL = (byte) Math.round(INDEL_ERROR_RATE * -10.0);

    /**
     * No indel likelihood (ref allele) used in the indel model to assess the confidence on the hom-ref call.
     */
    private static final double NO_INDEL_LIKELIHOOD = QualityUtils.qualToProbLog10(INDEL_QUAL);

    /**
     * Indel likelihood (alt. allele) used in the indel model to assess the confidence on the hom-ref call.
     */
    private static final double INDEL_LIKELIHOOD = QualityUtils.qualToErrorProbLog10(INDEL_QUAL);
    private static final int IDX_HOM_REF = 0;

    /**
     * Options related to posterior probability calcs
     */
    private static final boolean useInputSamplesAlleleCounts = false;  //by definition ref-conf will be single-sample; inputs should get ignored but let's be explicit
    private static final boolean useMLEAC = true;
    private static final boolean ignoreInputSamplesForMissingVariants = true;
    private static final boolean useFlatPriorsForIndels = false;


    /**
     * Create a new ReferenceConfidenceModel
     *
     * @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     */
    public ReferenceConfidenceModel(final SampleList samples,
                                    final SAMFileHeader header,
                                    final int indelInformativeDepthIndelSize,
                                    final int numRefForPrior) {
        Utils.nonNull(samples, "samples cannot be null");
        Utils.validateArg( samples.numberOfSamples() > 0, "samples cannot be empty");
        Utils.nonNull(header, "header cannot be empty");
        //TODO: code and comment disagree -- which is right?
        Utils.validateArg( indelInformativeDepthIndelSize >= 0, () -> "indelInformativeDepthIndelSize must be >= 1 but got " + indelInformativeDepthIndelSize);

        this.samples = samples;
        this.indelInformativeDepthIndelSize = indelInformativeDepthIndelSize;
        this.numRefSamplesForPrior = numRefForPrior;
        this.options = new PosteriorProbabilitiesUtils.PosteriorProbabilitiesOptions(HomoSapiensConstants.SNP_HETEROZYGOSITY,
                HomoSapiensConstants.INDEL_HETEROZYGOSITY, useInputSamplesAlleleCounts, useMLEAC, ignoreInputSamplesForMissingVariants,
                useFlatPriorsForIndels);
    }

    /**
     * Get the VCF header lines to include when emitting reference confidence values via {@link #calculateRefConfidence}.
     * @return a non-null set of VCFHeaderLines
     */
    public Set<VCFHeaderLine> getVCFHeaderLines() {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();
        headerLines.add(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME, "Represents any possible alternative allele at this location"));
        return headerLines;
    }

    public List<VariantContext> calculateRefConfidence(final Haplotype refHaplotype,
                                                       final Collection<Haplotype> calledHaplotypes,
                                                       final SimpleInterval paddedReferenceLoc,
                                                       final AssemblyRegion activeRegion,
                                                       final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                       final PloidyModel ploidyModel,
                                                       final List<VariantContext> variantCalls) {
        return calculateRefConfidence(refHaplotype, calledHaplotypes, paddedReferenceLoc, activeRegion, readLikelihoods,
                ploidyModel, variantCalls, false, Collections.emptyList());
    }

    /**
     * Calculate the reference confidence for a single sample given the its read data
     *
     * Returns a list of variant contexts, one for each position in the {@code activeRegion.getLoc()}, each containing
     * detailed information about the certainty that the sample is hom-ref for each base in the region.
     *
     *
     *
     * @param refHaplotype the reference haplotype, used to get the reference bases across activeRegion.getLoc()
     * @param calledHaplotypes a list of haplotypes that segregate in this region, for realignment of the reads in the
     *                         readLikelihoods, corresponding to each reads best haplotype.  Must contain the refHaplotype.
     * @param paddedReferenceLoc the location of refHaplotype (which might be larger than activeRegion.getLoc())
     * @param activeRegion the active region we want to get the reference confidence over
     * @param readLikelihoods a map from a single sample to its PerReadAlleleLikelihoodMap for each haplotype in calledHaplotypes
     * @param ploidyModel indicate the ploidy of each sample in {@code stratifiedReadMap}.
     * @param variantCalls calls made in this region.  The return result will contain any variant call in this list in the
     *                     correct order by genomic position, and any variant in this list will stop us emitting a ref confidence
     *                     under any position it covers (for snps and insertions that is 1 bp, but for deletions its the entire ref span)
     * @return an ordered list of variant contexts that spans activeRegion.getLoc() and includes both reference confidence
     *         contexts as well as calls from variantCalls if any were provided
     */
    public List<VariantContext> calculateRefConfidence(final Haplotype refHaplotype,
                                                       final Collection<Haplotype> calledHaplotypes,
                                                       final SimpleInterval paddedReferenceLoc,
                                                       final AssemblyRegion activeRegion,
                                                       final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                       final PloidyModel ploidyModel,
                                                       final List<VariantContext> variantCalls,
                                                       final boolean applyPriors,
                                                       final List<VariantContext> VCpriors) {
        Utils.nonNull(refHaplotype, "refHaplotype cannot be null");
        Utils.nonNull(calledHaplotypes, "calledHaplotypes cannot be null");
        Utils.validateArg(calledHaplotypes.contains(refHaplotype), "calledHaplotypes must contain the refHaplotype");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");
        Utils.nonNull(activeRegion, "activeRegion cannot be null");
        Utils.nonNull(readLikelihoods, "readLikelihoods cannot be null");
        Utils.validateArg(readLikelihoods.numberOfSamples() == 1, () -> "readLikelihoods must contain exactly one sample but it contained " + readLikelihoods.numberOfSamples());
        Utils.validateArg( refHaplotype.length() == activeRegion.getExtendedSpan().size(), () -> "refHaplotype " + refHaplotype.length() + " and activeRegion location size " + activeRegion.getSpan().size() + " are different");
        Utils.nonNull(ploidyModel, "the ploidy model cannot be null");
        final int ploidy = ploidyModel.samplePloidy(0); // the first sample = the only sample in reference-confidence mode.

        final SimpleInterval refSpan = activeRegion.getSpan();
        final List<ReadPileup> refPileups = AssemblyBasedCallerUtils.getPileupsOverReference(activeRegion.getHeader(), refSpan, readLikelihoods, samples);
        final byte[] ref = refHaplotype.getBases();
        final List<VariantContext> results = new ArrayList<>(refSpan.size());
        final String sampleName = readLikelihoods.getSample(0);

        final int globalRefOffset = refSpan.getStart() - activeRegion.getExtendedSpan().getStart();
        // Note, we use an indexed for-loop here because this method has a large impact on the profile of HaplotypeCaller runtime in GVCF mode
        final int refPileupsSize = refPileups.size();
        for (int i = 0; i < refPileupsSize; i++) {
            final ReadPileup pileup = refPileups.get(i);
            final Locatable curPos = pileup.getLocation();
            final int offset = curPos.getStart() - refSpan.getStart();

            final VariantContext overlappingSite = GATKVariantContextUtils.getOverlappingVariantContext(curPos, variantCalls);
            final List<VariantContext> currentPriors = VCpriors.isEmpty() ? Collections.emptyList() : getMatchingPriors(curPos, overlappingSite, VCpriors);
            if (overlappingSite != null && overlappingSite.getStart() == curPos.getStart()) {
                if (applyPriors) {
                    results.add(PosteriorProbabilitiesUtils.calculatePosteriorProbs(overlappingSite, currentPriors,
                            numRefSamplesForPrior, options));
                } else {
                    results.add(overlappingSite);
                }
            } else {
                // otherwise emit a reference confidence variant context
                results.add(makeReferenceConfidenceVariantContext(ploidy, ref, sampleName, globalRefOffset, pileup, curPos, offset, applyPriors, currentPriors));
            }
        }

        // Ensuring that we remove any indel informativeness data we may have attached to the underlying reads for caching purposes
        // This is important as if multiple reference blocks are computed for a low complexity active region some reads may incorrectly
        // be using caching values computed for a different reference block.
        if (USE_CACHED_READ_INDEL_INFORMATIVENESS_VALUES) {
            readLikelihoods.sampleEvidence(0).forEach(r -> r.clearTransientAttribute(INDEL_INFORMATIVE_BASES_CACHE_ATTRIBUTE_NAME));
        }

        return results;
    }


   public VariantContext makeReferenceConfidenceVariantContext(final int ploidy,
                                                                 final byte[] ref,
                                                                 final String sampleName,
                                                                 final int globalRefOffset,
                                                                 final ReadPileup pileup,
                                                                 final Locatable curPos,
                                                                 final int offset,
                                                                 final boolean applyPriors,
                                                                 final List<VariantContext> VCpriors) {
        // Assume infinite population on a single sample.
        final int refOffset = offset + globalRefOffset;
        final byte refBase = ref[refOffset];
        final ReferenceConfidenceResult homRefCalc = calcGenotypeLikelihoodsOfRefVsAny(ploidy, pileup, refBase, BASE_QUAL_THRESHOLD, null, true);

        final Allele refAllele = Allele.create(refBase, true);
        final List<Allele> refSiteAlleles = Arrays.asList(refAllele, Allele.NON_REF_ALLELE);
        final VariantContextBuilder vcb = new VariantContextBuilder("HC", curPos.getContig(), curPos.getStart(), curPos.getStart(), refSiteAlleles);
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, GATKVariantContextUtils.homozygousAlleleList(refAllele, ploidy));
        gb.AD(homRefCalc.getAD());
        gb.DP(homRefCalc.getDP());

        doIndelRefConfCalc(ploidy, ref, pileup, refOffset, homRefCalc);

       addGenotypeData(homRefCalc, gb);
        if(!applyPriors) {
            return vcb.genotypes(gb.make()).make();
        }
        else {
            return PosteriorProbabilitiesUtils.calculatePosteriorProbs(vcb.genotypes(gb.make()).make(), VCpriors, numRefSamplesForPrior, options);
            //TODO FIXME: after new-qual refactoring, these should be static calls to AF calculator
        }
    }

    public void doIndelRefConfCalc(final int ploidy, final byte[] ref, final ReadPileup pileup, final int refOffset, final ReferenceConfidenceResult refResult) {
        final RefVsAnyResult homRefCalc = (RefVsAnyResult)refResult;
        // genotype likelihood calculation
        final GenotypeLikelihoods snpGLs = GenotypeLikelihoods.fromLog10Likelihoods(homRefCalc.getGenotypeLikelihoodsCappedByHomRefLikelihood());
        final int nIndelInformativeReads = calcNReadsWithNoPlausibleIndelsReads(pileup, refOffset, ref, indelInformativeDepthIndelSize);
        final GenotypeLikelihoods indelGLs = getIndelPLs(ploidy,nIndelInformativeReads);

        // now that we have the SNP and indel GLs, we take the one with the least confidence,
        // as this is the most conservative estimate of our certainty that we are hom-ref.
        // For example, if the SNP PLs are 0,10,100 and the indel PLs are 0,100,1000
        // we are very certain that there's no indel here, but the SNP confidence imply that we are
        // far less confident that the ref base is actually the only thing here.  So we take 0,10,100
        // as our GLs for the site.
        final GenotypeLikelihoods leastConfidenceGLs = getGLwithWorstGQ(indelGLs, snpGLs);

        homRefCalc.finalPhredScaledGenotypeLikelihoods = leastConfidenceGLs.getAsPLs();
    }

    public void addGenotypeData(final ReferenceConfidenceResult result, final GenotypeBuilder gb) {
        final int[] pls = ((RefVsAnyResult)result).finalPhredScaledGenotypeLikelihoods;
        gb.PL(pls);
        gb.GQ(GATKVariantContextUtils.calculateGQFromPLs(pls));
    }

    /**
     * Get the GenotypeLikelihoods with the least strong corresponding GQ value
     * @param gl1 first to consider (cannot be null)
     * @param gl2 second to consider (cannot be null)
     * @return gl1 or gl2, whichever has the worst GQ
     */
    @VisibleForTesting
    GenotypeLikelihoods getGLwithWorstGQ(final GenotypeLikelihoods gl1, final GenotypeLikelihoods gl2) {
        if (getGQForHomRef(gl1) > getGQForHomRef(gl2)) {
            return gl1;
        } else {
            return gl2;
        }
    }

    private double getGQForHomRef(final GenotypeLikelihoods gls){
        return GenotypeLikelihoods.getGQLog10FromLikelihoods(IDX_HOM_REF, gls.getAsVector());
    }

    /**
     * Get indel PLs corresponding to seeing N nIndelInformativeReads at this site
     *
     * @param nInformativeReads the number of reads that inform us about being ref without an indel at this site
     * @param ploidy the requested ploidy.
     * @return non-null GenotypeLikelihoods given N
     */
    @VisibleForTesting
    GenotypeLikelihoods getIndelPLs(final int ploidy, final int nInformativeReads) {
        return indelPLCache(ploidy, nInformativeReads > MAX_N_INDEL_INFORMATIVE_READS ? MAX_N_INDEL_INFORMATIVE_READS : nInformativeReads);
    }

    private GenotypeLikelihoods indelPLCache(final int ploidy, final int nInformativeReads) {
        return initializeIndelPLCache(ploidy)[nInformativeReads];
    }

    private GenotypeLikelihoods[] initializeIndelPLCache(final int ploidy) {

        if (indelPLCache.length <= ploidy) {
            indelPLCache = Arrays.copyOf(indelPLCache, ploidy << 1);
        }

        if (indelPLCache[ploidy] != null) {
            return indelPLCache[ploidy];
        }

        final double denominator =  - MathUtils.log10(ploidy);
        final GenotypeLikelihoods[] result = new GenotypeLikelihoods[MAX_N_INDEL_INFORMATIVE_READS + 1];

        //Note: an array of zeros is the right answer for result[0].
        result[0] = GenotypeLikelihoods.fromLog10Likelihoods(new double[ploidy + 1]);
        for( int nInformativeReads = 1; nInformativeReads <= MAX_N_INDEL_INFORMATIVE_READS; nInformativeReads++ ) {
            final double[] PLs = new double[ploidy + 1];
            PLs[0] = nInformativeReads * NO_INDEL_LIKELIHOOD;
            for (int altCount = 1; altCount <= ploidy; altCount++) {
                final double refLikelihoodAccum = NO_INDEL_LIKELIHOOD + MathUtils.log10(ploidy - altCount);
                final double altLikelihoodAccum = INDEL_LIKELIHOOD + MathUtils.log10(altCount);
                PLs[altCount] = nInformativeReads * (MathUtils.approximateLog10SumLog10(refLikelihoodAccum ,altLikelihoodAccum) + denominator);
            }
            result[nInformativeReads] = GenotypeLikelihoods.fromLog10Likelihoods(PLs);
        }
        indelPLCache[ploidy] = result;
        return result;
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
    public ReferenceConfidenceResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
                                                        final ReadPileup pileup,
                                                        final byte refBase,
                                                        final byte minBaseQual,
                                                        final MathUtils.RunningAverage hqSoftClips,
                                                            final boolean readsWereRealigned) {

        final int likelihoodCount = ploidy + 1;
        final double log10Ploidy = MathUtils.log10(ploidy);

        final RefVsAnyResult result = new RefVsAnyResult(likelihoodCount);
        int readCount = 0;
        for (final PileupElement p : pileup) {
            final byte qual = p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual();
            if (!p.isDeletion() && qual <= minBaseQual) {
                continue;
            }
            readCount++;
            applyPileupElementRefVsNonRefLikelihoodAndCount(refBase, likelihoodCount, log10Ploidy, result, p, qual, hqSoftClips, readsWereRealigned);
        }
        final double denominator = readCount * log10Ploidy;
        for (int i = 0; i < likelihoodCount; i++) {
            result.genotypeLikelihoods[i] -= denominator;
        }
        return result;
    }

    private void applyPileupElementRefVsNonRefLikelihoodAndCount(final byte refBase, final int likelihoodCount, final double log10Ploidy, final RefVsAnyResult result, final PileupElement element, final byte qual, final MathUtils.RunningAverage hqSoftClips, final boolean readsWereRealigned) {
        final boolean isAlt = readsWereRealigned ? isAltAfterAssembly(element, refBase) : isAltBeforeAssembly(element, refBase);
        final double referenceLikelihood;
        final double nonRefLikelihood;
        if (isAlt) {
            nonRefLikelihood = QualityUtils.qualToProbLog10(qual);
            referenceLikelihood = QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG10_ONE_THIRD;
            result.nonRefDepth++;
        } else {
            referenceLikelihood = QualityUtils.qualToProbLog10(qual);
            nonRefLikelihood = QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG10_ONE_THIRD;
            result.refDepth++;
        }
        // Homozygous likelihoods don't need the logSum trick.
        result.genotypeLikelihoods[0] += referenceLikelihood + log10Ploidy;
        result.genotypeLikelihoods[likelihoodCount - 1] += nonRefLikelihood + log10Ploidy;
        // Heterozygous likelihoods need the logSum trick:
        for (int i = 1, j = likelihoodCount - 2; i < likelihoodCount - 1; i++, j--) {
            result.genotypeLikelihoods[i] +=
                    MathUtils.approximateLog10SumLog10(
                            referenceLikelihood + MathUtils.log10(j),
                            nonRefLikelihood + MathUtils.log10(i));
        }
        if (isAlt && hqSoftClips != null && element.isNextToSoftClip()) {
            hqSoftClips.add(AlignmentUtils.calcNumHighQualitySoftClips(element.getRead(), HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD));
        }
    }

    protected static boolean isAltBeforeAssembly(final PileupElement element, final byte refBase){
        return element.getBase() != refBase || element.isDeletion() || element.isBeforeDeletionStart()
                || element.isAfterDeletionEnd() || element.isBeforeInsertion() || element.isAfterInsertion() || element.isNextToSoftClip();
    }

    protected static boolean isAltAfterAssembly(final PileupElement element, final byte refBase){
        return element.getBase() != refBase || element.isDeletion(); //we shouldn't have soft clips after assembly
    }


    /**
     * Note that we don't have to match alleles because the PosteriorProbabilitesUtils will take care of that
     * @param curPos position of interest for genotyping
     * @param call (may be null)
     * @param priorList priors within the current ActiveRegion
     * @return prior VCs representing the same variant position as call
     */
    private List<VariantContext> getMatchingPriors(final Locatable curPos, final VariantContext call, final List<VariantContext> priorList) {
        final int position = call != null ? call.getStart() : curPos.getStart();
        final List<VariantContext> matchedPriors = new ArrayList<>(priorList.size());
        // NOTE: a for loop is used here because this method ends up being called per-pileup, per-read and using a loop instead of streaming saves runtime
        final int priorsListSize = priorList.size();
        for (int i = 0; i < priorsListSize; i++) {
            if (position == priorList.get(i).getStart()) {
                matchedPriors.add(priorList.get(i));
            }
        }
        return matchedPriors;
    }

    /**
     * Compute the sum of mismatching base qualities for readBases aligned to refBases at readStart / refStart
     * assuming no insertions or deletions in the read w.r.t. the reference
     *
     * @param readBases non-null bases of the read
     * @param readQuals non-null quals of the read
     * @param readStart the starting position of the read (i.e., that aligns it to a position in the reference)
     * @param refBases the reference bases
     * @param refStart the offset into refBases that aligns to the readStart position in readBases
     * @return an array containing the sum of quality scores for readBases that mismatch following this base and their corresponding ref base for each read base in readBases
     */
    private static int[] calculateBaselineMMQualities(final byte[] readBases,
                                final byte[] readQuals,
                                final int readStart,
                                final byte[] refBases,
                                final int refStart) {
        final int n = Math.min(readBases.length - readStart, refBases.length - refStart);
        int[] results = new int[n];
        int sum = 0;

        // Note that we start this loop at the end based on the principle that in order to calculate the number of mismatches remaining
        // between the read and the reference after the nth base, one can simply first calculate the remaining mismatches for the n + 1th
        // base first and so on.
        for ( int i = n - 1; i >= 0; i-- ) {
            final byte readBase = readBases[readStart + i];
            final byte refBase  = refBases[refStart + i];
            if (isMismatchAndNotAnAlignmentGap(readBase, refBase)) {
                sum += readQuals[readStart + i];
            }
            results[i] = sum;
        }

        return results;
    }

    /**
     * Compute whether a read is informative to eliminate an indel of size <= maxIndelSize segregating at readStart/refStart
     *
     * For each base this method determines if there are any plausible indels of size <= maxIndelSize that start at that
     * base. The method returns true if no indels were found that align as well or better than the rest of this read
     * compared to the reference.
     *
     * In the computation of this function for a given base, it also computes the value for every readOffset to the
     * end of the read as well. These results are cached in a bitset in the transient attributes for the read. A 1 in
     * the bitset means that this method would return true for that particular readOffset/refOffset combination. Note
     * that if a bitset is found in on the read already, this method defaults to returning the cached value over
     * computing the plausible indels again.
     *
     * Positions <= maxIndelSize from the end of the provided read/ref always return false.
     *
     * ***WARNING: the caching code makes the assumption that this function will be called over reference bases in ascending order. The
     *       results are undefined and will likely be wrong if used in any other way. If calling this method out of order, set
     *       useCachedResults to false
     *
     * @param read the read
     * @param readStart the 0-based index with respect to @{param}refBases where the read starts (this is the "IGV View" offset for the read)
     * @param refBases the reference bases
     * @param refStart the 0-based offset into refBases that aligns to the readStart position in readBases
     * @param maxIndelSize the max indel size to consider for the read to be informative
     * @param useCachedResults if false, ignore cached results for informative indel sizes (useful for debugging)
     * @return true if read can eliminate the possibility that there's an indel of size <= maxIndelSize segregating at refStart
     */
    private static boolean readHasNoPlausibleIdealsOfSize(final GATKRead read,
                                                          final int readStart,
                                                          final byte[] refBases,
                                                          final int refStart,
                                                          final int maxIndelSize,
                                                          final boolean useCachedResults) {
        BitSet cachedResult = (BitSet) read.getTransientAttribute(INDEL_INFORMATIVE_BASES_CACHE_ATTRIBUTE_NAME);
        if (cachedResult == null || !useCachedResults) {
            Utils.validate(readStart >= 0, "readStart must >= 0");
            Utils.validate(refStart >= 0, "refStart must >= 0");
            BitSet informativeBases = new BitSet(read.getLength());

            // Check that we aren't so close to the end of the end of the read that we don't have to compute anything more
            if ( !(read.getLength() - readStart < maxIndelSize) && !(refBases.length - refStart < maxIndelSize) ) {
                //TODO this should be removed, see https://github.com/broadinstitute/gatk/issues/5646 to track its progress
                final int secondaryReadBreakPosition = read.getLength() - maxIndelSize;

                // We are safe to use the faster no-copy versions of getBases and getBaseQualities here,
                // since we're not modifying the returned arrays in any way. This makes a small difference
                // in the HaplotypeCaller profile, since this method is a major hotspot.
                final Pair<byte[], byte[]> readBasesAndBaseQualities = AlignmentUtils.getBasesAndBaseQualitiesAlignedOneToOne(read);  //calls getBasesNoCopy if CIGAR is all match
                final byte[] readBases = readBasesAndBaseQualities.getLeft();
                final byte[] readQualities = readBasesAndBaseQualities.getRight();

                // Need to check for closeness to the end of the read again as the array size may be different than read.Len() due to deletions in the cigar
                if (readBases.length - readStart > maxIndelSize) {

                    // Compute where the end of marking would have been given the above two break conditions so we can stop marking there for our cached results
                    final int lastReadBaseToMarkAsIndelRelevant;
                    final boolean referenceWasShorter;
                    if (readBases.length < refBases.length - refStart + readStart + 1) {
                        // If the read ends first, then we don't mark the last maxIndelSize bases from it as relevant
                        lastReadBaseToMarkAsIndelRelevant = readBases.length - maxIndelSize;
                        referenceWasShorter = false;
                    } else {
                        // If the reference ends first, then we don't mark the last maxIndelSize bases from it as relevant
                        lastReadBaseToMarkAsIndelRelevant = refBases.length - refStart + readStart - maxIndelSize + 1;
                        referenceWasShorter = true;
                    }


                    // Compute the absolute baseline sum against which to test
                    final int[] baselineMisMatchSums = calculateBaselineMMQualities(readBases, readQualities, readStart, refBases, refStart);

                    // consider each indel size up to max in term, checking if an indel that deletes either the ref bases (deletion)
                    // or read bases (insertion) would fit as well as the origin baseline sum of mismatching quality scores. These scores
                    // are computed starting from the last base in the read/reference that would be offset by the indel and compared against
                    // the mismatch cost for the same base of the reference. Once the sum of mismatch qualities counting from the back for
                    // one indel size exceeds the global indel mismatch cost, the code stops as it will never find a better mismatch value.
                    for (int indelSize = 1; indelSize <= maxIndelSize; indelSize++) {
                        // Computing mismatches corresponding to a deletion
                        traverseEndOfReadForIndelMismatches(informativeBases,
                                readStart,
                                readBases,
                                readQualities,
                                lastReadBaseToMarkAsIndelRelevant,
                                secondaryReadBreakPosition,
                                refStart,
                                refBases,
                                baselineMisMatchSums,
                                indelSize,
                                false);

                        // Computing mismatches corresponding to an insertion
                        traverseEndOfReadForIndelMismatches(informativeBases,
                                readStart,
                                readBases,
                                readQualities,
                                lastReadBaseToMarkAsIndelRelevant,
                                secondaryReadBreakPosition,
                                refStart,
                                refBases,
                                baselineMisMatchSums,
                                indelSize,
                                true);
                    }


                    // Flip the bases at the front of the read (the ones not within maxIndelSize of the end as those are never informative)
                    // These must be flipped because thus far we have marked reads for which there were plausible indels with a true value in
                    // the bitset. This method returns false for cases where we have discovered plausible indels so we must flip them. This
                    // is done in part to preserve a sensible default behavior for bases not considered by this approach.
                    if ( lastReadBaseToMarkAsIndelRelevant <= secondaryReadBreakPosition) {
                        informativeBases.flip(0, lastReadBaseToMarkAsIndelRelevant);
                        // Resolve the fact that the old approach would always mark the last base examined as being indel uninformative when the reference
                        // ends first despite it corresponding to a comparison of zero bases against the read
                        if (referenceWasShorter) {
                            informativeBases.set(lastReadBaseToMarkAsIndelRelevant - 1, false);
                        }
                    } else {
                        informativeBases.flip(0, secondaryReadBreakPosition + 1);
                    }

                }
            }
            cachedResult = informativeBases;
            read.setTransientAttribute(INDEL_INFORMATIVE_BASES_CACHE_ATTRIBUTE_NAME, informativeBases);
        }
        return cachedResult.get(readStart);
    }

    /**
     * Helper method responsible for read-end traversal. This method will handle both insertions and deletions,
     * indicated by setting the insertion parameter.
     *
     * Given the array of sums baselineMMSums, this method will start from the back of the read and reference and sum the
     * quality score of all mismatching bases to the reference. If the score is equal to or lower than the baseline sum and
     * the base being examined is before lastReadBaseToMarkAsIndelRelevant, then this method will store a true into the informativeBases
     * Bitset for that particular read. If at any point the sum for a given size insertion/deletion exceeds the global cost
     * of all aligned mismatches to the reference with no indels added (the first position in baselineMMSums) the process will
     * end prematurely so as to avoid comparing any additional bases beyond what is necessary.
     *
     * It is expected that only the bases between readStart and lastReadBaseToMarkAsIndelRelevant in the bitset will be set to true
     * by this method if they are ambiguous about an indel of the given size. We then flip these values later in the process because
     * an ambiguous indel positions in the read actually return false in readHasNoPlausibleIdealsOfSize.
     *
     * NOTE: This method examines overhanging bases to the reference/read if they do not end at the same position.
     *       (eg. if the reference ends 20 bases after the read does and you are looking at a deletion of size 5, the first
     *       base compared will be the last base of the read and the 15th from last base on the reference)
     *
     * @param informativeBases ReadBases indexed bitset into which to store the results
     * @param readStart Offset of first comparison base into the read
     * @param readBases Read bases aligned to be indexed by reference base
     * @param readQuals Read qualities aligned to be indexed by reference base
     * @param lastReadBaseToMarkAsIndelRelevant Final base in the read that is valid to store results for based on closeness to the edge of the indel.
     *                                          This method may compare bases beyond this point but it will not mark them as being relevant in the output
     *                                          unless the read base being compared lies before this index.
     * @param secondaryReadBreakPosition Break position to compare based on the read.Length() ending position compared to the readBases.length
     * @param refStart Starting base in the reference array to consider
     * @param refBases Array of reference bases to compare
     * @param baselineMMSums Array of mismatch scores for each position on the read. (NOTE this array should not be mutated
     *                       by this method as it is shared between calls to this method)
     * @param indelSize size of offset between reference and read bases to consider
     * @param insertion whether to compute offsets for an insertion (otherwise treats the offset as a deletion)
     */
    private static void traverseEndOfReadForIndelMismatches(final BitSet informativeBases, final int readStart, final byte[] readBases, final byte[] readQuals, final int lastReadBaseToMarkAsIndelRelevant,  final int secondaryReadBreakPosition, final int refStart, final byte[] refBases,  final int[] baselineMMSums, final int indelSize, final boolean insertion) {
        final int globalMismatchCostForReadAlignedToReference = baselineMMSums[0];
        int baseQualitySum = 0;

        // Compute how many bases forward we should compare taking into account reference/read overhang
        final int insertionLength = !insertion ? 0 : indelSize;
        final int deletionLength = insertion ? 0 : indelSize;

        // Based on the offsets and the indelSize we are considering, how many bases until we fall off the end of the read/reference arrays?
        final int numberOfBasesToDirectlyCompare = Math.min(readBases.length - readStart - insertionLength,
                refBases.length - refStart - deletionLength);

        for (int readOffset = numberOfBasesToDirectlyCompare + insertionLength - 1,
             refOffset = numberOfBasesToDirectlyCompare + deletionLength - 1;
             readOffset >= 0 && refOffset >= 0;
             readOffset--, refOffset--) {

            // Calculate the real base offset for the read:
            final byte readBase = readBases[readStart + readOffset];
            final byte refBase = refBases[refStart + refOffset];
            if (isMismatchAndNotAnAlignmentGap(readBase, refBase)) {
                baseQualitySum += readQuals[readStart + readOffset];
                if (baseQualitySum > globalMismatchCostForReadAlignedToReference) { // abort early if we are over our global mismatch cost
                    break;
                }
            }
            // The hypothetical "readOffset" that corresponds to the comparison we are currently making
            int siteOfRealComparisonPoint = Math.min(readOffset, refOffset);

            // If it's a real character and the cost isn't greater than the non-indel cost, label it as uninformative
            if (readBases[readStart + siteOfRealComparisonPoint] != AlignmentUtils.GAP_CHARACTER &&
                    // Use less than here because lastReadBaseToMarkAsIndelRelevant is the exclusive site where we flip bases later on.
                    readStart + siteOfRealComparisonPoint < lastReadBaseToMarkAsIndelRelevant &&
                    // Resolving the edge case involving read.getLength() disagreeing with the realigned indel length
                    readStart + siteOfRealComparisonPoint <= secondaryReadBreakPosition &&
                    baselineMMSums[siteOfRealComparisonPoint] >= baseQualitySum) {
                informativeBases.set(readStart + siteOfRealComparisonPoint, true); // Label with true here because we flip these results later
            }
        }
    }

    // Are these two bases different (including IUPAC bases) and does the read not correspond to a deletion on the reference
    private static boolean isMismatchAndNotAnAlignmentGap(byte readBase, byte refBase) {
        return !Nucleotide.intersect(readBase, refBase) && (readBase != AlignmentUtils.GAP_CHARACTER);
    }

    /**
     * Calculate the number of reads that have no plausible indels based on alignment at pileup
     *
     * @param pileup a pileup
     * @param pileupOffsetIntoRef index along the reference corresponding to the pileup
     * @param ref the ref bases
     * @param maxIndelSize maximum indel size to consider in the as plausible calculation
     * @return an integer >= 0
     */
    @VisibleForTesting
    int calcNReadsWithNoPlausibleIndelsReads(final ReadPileup pileup, final int pileupOffsetIntoRef, final byte[] ref, final int maxIndelSize) {
        int nInformative = 0;
        for ( final PileupElement p : pileup ) {
            // doesn't count as evidence
            if ( p.isBeforeDeletionStart() || p.isBeforeInsertion() || p.isDeletion() ) {
                continue;
            }

            final int offset = getCigarModifiedOffset(p);

            if ( readHasNoPlausibleIdealsOfSize(p.getRead(), offset, ref, pileupOffsetIntoRef, maxIndelSize, USE_CACHED_READ_INDEL_INFORMATIVENESS_VALUES) ) {
                nInformative++;
                if( nInformative > MAX_N_INDEL_INFORMATIVE_READS ) {
                    return MAX_N_INDEL_INFORMATIVE_READS;
                }
            }
        }
        return nInformative;
    }

    /**
     * Calculate the index of the current pileup position against the reference-aligned read
     * This offset should be representative of the "IGV view" for the read where insertions are collapsed and deletions
     * are padded so that we can easily count the mismatches against the reference
     * @param p the PileupElement containing the offset as an index into the read base sequence
     * @return the new reference-aligned index/offset
     */
    @VisibleForTesting
    protected int getCigarModifiedOffset (final PileupElement p){
        final GATKRead read = p.getRead();
        int offset = (p.getCurrentCigarElement().getOperator().consumesReferenceBases() || p.getCurrentCigarElement().getOperator() == CigarOperator.S)? p.getOffsetInCurrentCigar() : 0;
        for (int i = 0; i < p.getCurrentCigarOffset(); i++) {
            final CigarElement elem = read.getCigarElement(i);
            if (elem.getOperator().consumesReferenceBases() || elem.getOperator() == CigarOperator.S) {
                offset += elem.getLength();
            }
        }
        return offset;
    }

    /**
     * Create a reference haplotype for an active region
     *
     * @param activeRegion the active region
     * @param refBases the ref bases
     * @param paddedReferenceLoc the location spanning of the refBases -- can be longer than activeRegion.getLocation()
     * @return a reference haplotype
     */
    public static Haplotype createReferenceHaplotype(final AssemblyRegion activeRegion, final byte[] refBases, final SimpleInterval paddedReferenceLoc) {
        Utils.nonNull(activeRegion, "null region");
        Utils.nonNull(refBases, "null refBases");
        Utils.nonNull(paddedReferenceLoc, "null paddedReferenceLoc");

        final int alignmentStart = activeRegion.getExtendedSpan().getStart() - paddedReferenceLoc.getStart();
        if ( alignmentStart < 0 ) {
            throw new IllegalStateException("Bad alignment start in createReferenceHaplotype " + alignmentStart);
        }
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        refHaplotype.setAlignmentStartHapwrtRef(alignmentStart);
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);
        return refHaplotype;
    }
}
