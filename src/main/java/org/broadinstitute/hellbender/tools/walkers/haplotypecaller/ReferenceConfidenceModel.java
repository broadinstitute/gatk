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
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel;
import org.broadinstitute.hellbender.tools.walkers.variantutils.PosteriorProbabilitiesUtils;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Code for estimating the reference confidence
 *
 * This code can estimate the probability that the data for a single sample is consistent with a
 * well-determined REF/REF diploid genotype.
 *
 */
public final class ReferenceConfidenceModel {

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
     * Holds information about a genotype call of a single sample reference vs. any non-ref event
     *
     * IMPORTANT PERFORMANCE NOTE!!! Allowing direct field access (within this class only) speeds up
     * the HaplotypeCaller by ~10% vs. accessing the fields indirectly via setters, as seen in a profiler.
     */
    @VisibleForTesting
    static final class RefVsAnyResult {
        /**
         * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
         *
         * Fields are visible because direct field access for this particular class has a major performance
         * impact on the HaplotypeCaller, as noted above, and the class itself is nested within
         * ReferenceConfidenceModel anyway.
         */
        final double[] genotypeLikelihoods;

        int refDepth = 0;
        int nonRefDepth = 0;

        /**
         * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
         * @param likelihoodCapacity the required capacity of the likelihood array, should match the possible number of
         *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
         *                           model non-sense genotype count if applies.
         * @throws IllegalArgumentException if {@code likelihoodCapacity} is negative.
         */
        public RefVsAnyResult(final int likelihoodCapacity) {
            ParamUtils.isPositiveOrZero(likelihoodCapacity, "likelihood capacity is negative");
            genotypeLikelihoods = new double[likelihoodCapacity];
        }
        
        /**
         * @return Get the DP (sum of AD values)
         */
        int getDP() {
            return refDepth + nonRefDepth;
        }

        /**
         * Return the AD fields. Returns a newly allocated array every time.
         */
        int[] getAD() {
            return new int[]{refDepth, nonRefDepth};
        }

        /**
         * Returns (a copy of) the array of genotype likelihoods
         * Caps the het and hom var likelihood values by the hom ref likelihood.
         * The capping is done on the fly.
         */
        double[] getGenotypeLikelihoodsCappedByHomRefLikelihood() {
            return DoubleStream.of(genotypeLikelihoods).map(d -> Math.min(d, genotypeLikelihoods[0])).toArray();
        }
    }

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
                                                       final ReadLikelihoods<Haplotype> readLikelihoods,
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
                                                       final ReadLikelihoods<Haplotype> readLikelihoods,
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
        final List<ReadPileup> refPileups = getPileupsOverReference(refHaplotype, calledHaplotypes, paddedReferenceLoc, activeRegion, refSpan, readLikelihoods);
        final byte[] ref = refHaplotype.getBases();
        final List<VariantContext> results = new ArrayList<>(refSpan.size());
        final String sampleName = readLikelihoods.getSample(0);

        final int globalRefOffset = refSpan.getStart() - activeRegion.getExtendedSpan().getStart();
        for ( final ReadPileup pileup : refPileups ) {
            final Locatable curPos = pileup.getLocation();
            final int offset = curPos.getStart() - refSpan.getStart();

            final VariantContext overlappingSite = getOverlappingVariantContext(curPos, variantCalls);
            final List<VariantContext> currentPriors = getMatchingPriors(curPos, overlappingSite, VCpriors);
            if ( overlappingSite != null && overlappingSite.getStart() == curPos.getStart() ) {
                if (applyPriors) {
                    results.add(PosteriorProbabilitiesUtils.calculatePosteriorProbs(overlappingSite, currentPriors,
                            numRefSamplesForPrior, options));
                }
                else {
                    results.add(overlappingSite);
                }
            } else {
                // otherwise emit a reference confidence variant context
                results.add(makeReferenceConfidenceVariantContext(ploidy, ref, sampleName, globalRefOffset, pileup, curPos, offset, applyPriors, currentPriors));
            }
        }

        return results;
    }


   private VariantContext makeReferenceConfidenceVariantContext(final int ploidy,
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
        final RefVsAnyResult homRefCalc = calcGenotypeLikelihoodsOfRefVsAny(ploidy, pileup, refBase, BASE_QUAL_THRESHOLD, null, true);

        final Allele refAllele = Allele.create(refBase, true);
        final List<Allele> refSiteAlleles = Arrays.asList(refAllele, Allele.NON_REF_ALLELE);
        final VariantContextBuilder vcb = new VariantContextBuilder("HC", curPos.getContig(), curPos.getStart(), curPos.getStart(), refSiteAlleles);
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, GATKVariantContextUtils.homozygousAlleleList(refAllele, ploidy));
        gb.AD(homRefCalc.getAD());
        gb.DP(homRefCalc.getDP());

        // genotype likelihood calculation
        final GenotypeLikelihoods snpGLs = GenotypeLikelihoods.fromLog10Likelihoods(homRefCalc.getGenotypeLikelihoodsCappedByHomRefLikelihood());
        final int nIndelInformativeReads = calcNIndelInformativeReads(pileup, refOffset, ref, indelInformativeDepthIndelSize);
        final GenotypeLikelihoods indelGLs = getIndelPLs(ploidy,nIndelInformativeReads);

        // now that we have the SNP and indel GLs, we take the one with the least confidence,
        // as this is the most conservative estimate of our certainty that we are hom-ref.
        // For example, if the SNP PLs are 0,10,100 and the indel PLs are 0,100,1000
        // we are very certain that there's no indel here, but the SNP confidence imply that we are
        // far less confident that the ref base is actually the only thing here.  So we take 0,10,100
        // as our GLs for the site.
        final GenotypeLikelihoods leastConfidenceGLs = getGLwithWorstGQ(indelGLs, snpGLs);

        final int[] leastConfidenceGLsAsPLs = leastConfidenceGLs.getAsPLs();
        gb.GQ(GATKVariantContextUtils.calculateGQFromPLs(leastConfidenceGLsAsPLs));
        gb.PL(leastConfidenceGLsAsPLs);

        if(!applyPriors) {
            return vcb.genotypes(gb.make()).make();
        }
        else {
            return PosteriorProbabilitiesUtils.calculatePosteriorProbs(vcb.genotypes(gb.make()).make(), VCpriors, numRefSamplesForPrior, options);
            //TODO FIXME: after new-qual refactoring, these should be static calls to AF calculator
        }
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
    public RefVsAnyResult calcGenotypeLikelihoodsOfRefVsAny(final int ploidy,
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

    private boolean isAltBeforeAssembly(final PileupElement element, final byte refBase){
        return element.getBase() != refBase || element.isDeletion() || element.isBeforeDeletionStart()
                || element.isAfterDeletionEnd() || element.isBeforeInsertion() || element.isAfterInsertion() || element.isNextToSoftClip();
    }

    private boolean isAltAfterAssembly(final PileupElement element, final byte refBase){
        return element.getBase() != refBase || element.isDeletion(); //we shouldn't have soft clips after assembly
    }

    /**
     * Get a list of pileups that span the entire active region span, in order, one for each position
     */
    private List<ReadPileup> getPileupsOverReference(final Haplotype refHaplotype,
                                                           final Collection<Haplotype> calledHaplotypes,
                                                           final SimpleInterval paddedReferenceLoc,
                                                           final AssemblyRegion activeRegion,
                                                           final SimpleInterval activeRegionSpan,
                                                           final ReadLikelihoods<Haplotype> readLikelihoods) {
        Utils.validateArg(calledHaplotypes.contains(refHaplotype), "calledHaplotypes must contain the refHaplotype");
        Utils.validateArg(readLikelihoods.numberOfSamples() == 1, () -> "readLikelihoods must contain exactly one sample but it contained " + readLikelihoods.numberOfSamples());

        final List<GATKRead> reads = new ArrayList<>(readLikelihoods.sampleReads(0));
        reads.sort(new ReadCoordinateComparator(activeRegion.getHeader()));  //because we updated the reads based on the local realignments we have to re-sort or the pileups will be... unpredictable

        final LocusIteratorByState libs = new LocusIteratorByState(reads.iterator(), LocusIteratorByState.NO_DOWNSAMPLING,
                false, samples.asSetOfSamples(), activeRegion.getHeader(), true);

        final int startPos = activeRegionSpan.getStart();
        final List<ReadPileup> pileups = new ArrayList<>(activeRegionSpan.getEnd() - startPos);
        AlignmentContext next = libs.advanceToLocus(startPos, true);
        for ( int curPos = startPos; curPos <= activeRegionSpan.getEnd(); curPos++ ) {
            if ( next != null && next.getLocation().getStart() == curPos ) {
                pileups.add(next.getBasePileup());
                next = libs.hasNext() ? libs.next() : null;
            } else {
                // no data, so we create empty pileups
                pileups.add(new ReadPileup(new SimpleInterval(activeRegionSpan.getContig(), curPos, curPos)));
            }
        }

        return pileups;
    }

    /**
     * Return the rightmost variant context in maybeOverlapping that overlaps curPos
     *
     * @param curPos non-null genome loc
     * @param maybeOverlapping a collection of variant contexts that might overlap curPos
     * @return a VariantContext, or null if none overlaps
     */
    @VisibleForTesting
    VariantContext getOverlappingVariantContext(final Locatable curPos, final Collection<VariantContext> maybeOverlapping) {
        final SimpleInterval curPosSI = new SimpleInterval(curPos);
        VariantContext overlaps = null;
        for ( final VariantContext vc : maybeOverlapping ) {
            if ( curPosSI.overlaps(vc) ) {
                if ( overlaps == null || vc.getStart() > overlaps.getStart() ) {
                    overlaps = vc;
                }
            }
        }
        return overlaps;
    }

    /**
     * Note that we don't have to match alleles because the PosteriorProbabilitesUtils will take care of that
     * @param curPos position of interest for genotyping
     * @param call (may be null)
     * @param priorList priors within the current ActiveRegion
     * @return prior VCs representing the same variant position as call
     */
    List<VariantContext> getMatchingPriors(final Locatable curPos, final VariantContext call, final List<VariantContext> priorList) {
        final int position = call != null ? call.getStart() : curPos.getStart();
        return priorList.stream().filter(vc -> position == vc.getStart()).collect(Collectors.toList());
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
     * @param maxSum if the sum goes over this value, return immediately
     * @return the sum of quality scores for readBases that mismatch their corresponding ref bases
     */
    @VisibleForTesting
    int sumMismatchingQualities(final byte[] readBases,
                                final byte[] readQuals,
                                final int readStart,
                                final byte[] refBases,
                                final int refStart,
                                final int maxSum) {
        final int n = Math.min(readBases.length - readStart, refBases.length - refStart);
        int sum = 0;

        for ( int i = 0; i < n; i++ ) {
            final byte readBase = readBases[readStart + i];
            final byte refBase  = refBases[refStart + i];
            if ( !Nucleotide.intersect(readBase, refBase) && !(readBase == AlignmentUtils.GAP_CHARACTER)) {
                sum += readQuals[readStart + i];
                if ( sum > maxSum ) { // abort early
                    return sum;
                }
            }
        }

        return sum;
    }

    /**
     * Compute whether a read is informative to eliminate an indel of size <= maxIndelSize segregating at readStart/refStart
     *
     * @param read the read
     * @param readStart the index with respect to @{param}refBases where the read starts
     * @param refBases the reference bases
     * @param refStart the offset into refBases that aligns to the readStart position in readBases
     * @param maxIndelSize the max indel size to consider for the read to be informative
     * @return true if read can eliminate the possibility that there's an indel of size <= maxIndelSize segregating at refStart
     */
    @VisibleForTesting
    boolean isReadInformativeAboutIndelsOfSize(final GATKRead read,
                                             final int readStart,
                                             final byte[] refBases,
                                             final int refStart,
                                             final int maxIndelSize) {
        // fast exit when n bases left < maxIndelSize
        if( read.getLength() - readStart < maxIndelSize || refBases.length - refStart < maxIndelSize ) {
            return false;
        }

        // We are safe to use the faster no-copy versions of getBases and getBaseQualities here,
        // since we're not modifying the returned arrays in any way. This makes a small difference
        // in the HaplotypeCaller profile, since this method is a major hotspot.
        final byte[] readBases = AlignmentUtils.getBasesAlignedOneToOne(read);  //calls getBasesNoCopy if CIGAR is all match
        final byte[] readQuals = AlignmentUtils.getBaseQualsAlignedOneToOne(read);


        final int baselineMMSum = sumMismatchingQualities(readBases, readQuals, readStart, refBases, refStart, Integer.MAX_VALUE);

        // consider each indel size up to max in term, checking if an indel that deletes either the ref bases (deletion
        // or read bases (insertion) would fit as well as the origin baseline sum of mismatching quality scores
        for ( int indelSize = 1; indelSize <= maxIndelSize; indelSize++ ) {
            // check insertions:
            if (sumMismatchingQualities(readBases, readQuals, readStart + indelSize, refBases, refStart, baselineMMSum) <= baselineMMSum) {
                return false;
            }
            // check deletions:
            if (sumMismatchingQualities(readBases, readQuals, readStart, refBases, refStart + indelSize, baselineMMSum) <= baselineMMSum) {
                return false;
            }
        }

        return true;
    }

    /**
     * Calculate the number of indel informative reads at pileup
     *
     * @param pileup a pileup
     * @param pileupOffsetIntoRef index along the reference corresponding to the pileup
     * @param ref the ref bases
     * @param maxIndelSize maximum indel size to consider in the informativeness calculation
     * @return an integer >= 0
     */
    @VisibleForTesting
    int calcNIndelInformativeReads(final ReadPileup pileup, final int pileupOffsetIntoRef, final byte[] ref, final int maxIndelSize) {
        int nInformative = 0;
        for ( final PileupElement p : pileup ) {
            // doesn't count as evidence
            if ( p.isBeforeDeletionStart() || p.isBeforeInsertion() || p.isDeletion() ) {
                continue;
            }

            final int offset = getCigarModifiedOffset(p);

            if ( isReadInformativeAboutIndelsOfSize(p.getRead(), offset, ref, pileupOffsetIntoRef, maxIndelSize) ) {
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
        for (final CigarElement elem : read.getCigar().getCigarElements().subList(0, p.getCurrentCigarOffset())) {
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
