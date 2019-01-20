package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.primitives.Doubles;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by David Benjamin on 9/15/16.
 */
public class Mutect2FilteringEngine {
    public static final double MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT = 0.9;
    private static final int IMPUTED_NORMAL_BASE_QUALITY = 30;  // only used if normal base quality annotation fails somehow
    private static final double MIN_NORMAL_ARTIFACT_RATIO = 0.1;    // don't call normal artifact if allele fraction in normal is much smaller than allele fraction in tumor
    private M2FiltersArgumentCollection MTFAC;
    private final Map<String, Double> contaminationBySample;
    private final double somaticPriorProb;
    private final Set<String> normalSamples;
    final Map<String, OverlapDetector<MinorAlleleFractionRecord>> tumorSegments;
    public static final String FILTERING_STATUS_VCF_KEY = "filtering_status";

    public Mutect2FilteringEngine(final M2FiltersArgumentCollection MTFAC, final Set<String> normalSamples, final Map<String, Double> contaminationBySample) {
        this.MTFAC = MTFAC;
        this.contaminationBySample = contaminationBySample;
        this.normalSamples = normalSamples;
        somaticPriorProb = Math.pow(10, MTFAC.log10PriorProbOfSomaticEvent);
        tumorSegments = MTFAC.tumorSegmentationTables.stream()
                .map(MinorAlleleFractionRecord::readFromFile)
                .collect(Collectors.toMap(ImmutablePair::getLeft, p -> OverlapDetector.create(p.getRight())));
    }

    private void applyContaminationFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {

        final List<ImmutablePair<Integer, Double>> depthsAndPosteriors = new ArrayList<>();

        for (final Genotype tumorGenotype : vc.getGenotypes()) {
            if (normalSamples.contains(tumorGenotype.getSampleName())) {
                continue;
            }

            final double contamination = contaminationBySample.get(tumorGenotype.getSampleName());

            final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(tumorGenotype, VCFConstants.ALLELE_FREQUENCY_KEY,
                    () -> new double[] {1.0}, 1.0);
            final int maxFractionIndex = MathUtils.maxElementIndex(alleleFractions);
            final int[] ADs = tumorGenotype.getAD();
            final int altCount = ADs[maxFractionIndex + 1];   // AD is all alleles, while AF is alts only, hence the +1 offset
            final int depth = (int) MathUtils.sum(ADs);
            final double[] negativeLog10AlleleFrequencies = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc,
                    GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE, () -> new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}, Double.POSITIVE_INFINITY);
            final double alleleFrequency = MathUtils.applyToArray(negativeLog10AlleleFrequencies, x -> Math.pow(10,-x))[maxFractionIndex];

            final double somaticLikelihood = 1.0 / (depth + 1);

            final double singleContaminantLikelihood = 2 * alleleFrequency * (1 - alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination /2)
                    + MathUtils.square(alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination);
            final double manyContaminantLikelihood = MathUtils.binomialProbability(depth, altCount, contamination * alleleFrequency);
            final double contaminantLikelihood = Math.max(singleContaminantLikelihood, manyContaminantLikelihood);
            final double posteriorProbOfContamination = (1 - somaticPriorProb) * contaminantLikelihood / ((1 - somaticPriorProb) * contaminantLikelihood + somaticPriorProb * somaticLikelihood);

            depthsAndPosteriors.add(ImmutablePair.of(altCount, posteriorProbOfContamination));
        }

        double posteriorProbOfContamination = weightedMedianPosteriorProbability(depthsAndPosteriors);

        filterResult.addAttribute(GATKVCFConstants.CONTAMINATION_QUAL_ATTRIBUTE, QualityUtils.errorProbToQual(posteriorProbOfContamination));
        if (posteriorProbOfContamination > MTFAC.maxContaminationProbability) {
            filterResult.addFilter(GATKVCFConstants.CONTAMINATION_FILTER_NAME);
        }
    }

    // weighted median -- what's the lowest posterior probability that accounts for samples with half of the total alt depth
    private double weightedMedianPosteriorProbability(List<ImmutablePair<Integer, Double>> depthsAndPosteriors) {
        final int totalAltDepth = depthsAndPosteriors.stream().mapToInt(ImmutablePair::getLeft).sum();

        // sort from lowest to highest posterior probability of artifact
        depthsAndPosteriors.sort(Comparator.comparingDouble(p -> p.getRight()));

        int cumulativeAltCount = 0;

        for (final ImmutablePair<Integer, Double> pair : depthsAndPosteriors) {
            cumulativeAltCount += pair.getLeft();
            if (cumulativeAltCount * 2 >= totalAltDepth) {
                return pair.getRight();
            }
        }
        return 0;
    }

    private void applyTriallelicFilter(final VariantContext vc, final FilterResult filterResult) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY)) {
            final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
            final long numPassingAltAlleles = Arrays.stream(tumorLods).filter(x -> x > MTFAC.TUMOR_LOD_THRESHOLD).count();

            if (numPassingAltAlleles > MTFAC.numAltAllelesThreshold) {
                filterResult.addFilter(GATKVCFConstants.MULTIALLELIC_FILTER_NAME);
            }
        }
    }

    private void applySTRFilter(final VariantContext vc, final FilterResult filterResult) {
        // STR contractions, such as ACTACTACT -> ACTACT, are overwhelmingly false positives so we hard filter by default
        if (vc.isIndel()) {
            final int[] rpa = vc.getAttributeAsList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY).stream()
                    .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
            final String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");
            if (rpa == null || rpa.length < 2 || ru.length() == 0) {
                return;
            }
            final int referenceSTRBaseCount = ru.length() * rpa[0];
            final int numPCRSlips = rpa[0] - rpa[1];
            if (referenceSTRBaseCount >= MTFAC.minPcrSlippageBases && Math.abs(numPCRSlips) == 1) {
                // calculate the p-value that out of n reads we would have at least k slippage reads
                // if this p-value is small we keep the variant (reject the PCR slippage hypothesis)
                final int[] ADs = sumADsOverSamples(vc, true, false);

                final int depth = ADs == null ? 0 : (int) MathUtils.sum(ADs);
                final double oneSidedPValueOfSlippage = (ADs == null || ADs.length < 2) ? 1.0 :
                        new BinomialDistribution(null, depth, MTFAC.pcrSlippageRate).cumulativeProbability(ADs[1] - 1, depth);
                if (oneSidedPValueOfSlippage > MTFAC.pcrSlippagePValueThreshold) {
                    filterResult.addFilter(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME);
                }
            }
        }
    }

    private int[] sumADsOverSamples(final VariantContext vc, final boolean includeTumor, final boolean includeNormal) {
        final int[] ADs = new int[vc.getNAlleles()];
        vc.getGenotypes().stream()
                .filter(g -> includeTumor || normalSamples.contains(g.getSampleName()))
                .filter(g -> includeNormal || !normalSamples.contains(g.getSampleName()))
                .map(Genotype::getAD).forEach(ad -> new IndexRange(0, vc.getNAlleles()).forEach(n -> ADs[n] += ad[n]));
        return ADs;
    }

    private double[] weightedAverageOfTumorAFs(final VariantContext vc) {
        double totalWeight = 0.0;
        final double[] AFs = new double[vc.getNAlleles() - 1];
        for (final Genotype g : vc.getGenotypes()) {
            if (normalSamples.contains(g.getSampleName())) {
                continue;
            } else {
                final double weight = MathUtils.sum(g.getAD());
                totalWeight += weight;
                final double[] sampleAFs = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                        () -> new double[] {0.0}, 0.0);
                MathArrays.scaleInPlace(weight, sampleAFs);
                MathUtils.addToArrayInPlace(AFs, sampleAFs);
            }
        }
        MathArrays.scaleInPlace(1/totalWeight, AFs);
        return AFs;
    }

    private static void applyPanelOfNormalsFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        final boolean siteInPoN = vc.hasAttribute(GATKVCFConstants.IN_PON_VCF_ATTRIBUTE);
        if (siteInPoN) {
            filterResult.addFilter(GATKVCFConstants.PON_FILTER_NAME);
        }
    }

    private void applyBaseQualityFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (!vc.hasAttribute(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY)) {
            return;
        }

        final List<Integer> baseQualityByAllele = vc.getAttributeAsIntList(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY, 0);
        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        if (baseQualityByAllele.get(indexOfMaxTumorLod + 1) < MTFAC.minMedianBaseQuality) {
            filterResult.addFilter(GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME);
        }
    }

    private void applyMappingQualityFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (!vc.hasAttribute(MappingQuality.KEY)) {
            return;
        }

        final List<Integer> indelLengths = vc.getIndelLengths();
        final int indelLength = indelLengths == null ? 0 : indelLengths.stream().mapToInt(Math::abs).max().orElseGet(() -> 0);
        final List<Integer> mappingQualityByAllele = vc.getAttributeAsIntList(MappingQuality.KEY, 0);

        // we use the mapping quality annotation of the alt allele in most cases, but for long indels we use the reference
        // annotation.  We have to do this because the indel, even if it maps uniquely, gets a poor mapping quality
        // by virtue of its mismatch.  The reference mapping quality is a decent proxy for the region's mappability.
        if (mappingQualityByAllele.get(indelLength < MTFAC.longIndelLength ? 1 : 0) < MTFAC.minMedianMappingQuality) {
            filterResult.addFilter(GATKVCFConstants.MEDIAN_MAPPING_QUALITY_FILTER_NAME);
        }
    }

    private void applyMedianFragmentLengthDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (!vc.hasAttribute(FragmentLength.KEY)) {
            return;
        }

        final List<Integer> fragmentLengthByAllele = vc.getAttributeAsIntList(FragmentLength.KEY, 0);

        if (Math.abs(fragmentLengthByAllele.get(1) - fragmentLengthByAllele.get(0)) > MTFAC.maxMedianFragmentLengthDifference) {
            filterResult.addFilter(GATKVCFConstants.MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME);
        }
    }

    private void applyReadPositionFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (!vc.hasAttribute(ReadPosition.KEY)) {
            return;
        }

        final List<Integer> readPositionByAllele = vc.getAttributeAsIntList(ReadPosition.KEY, 0);

        // a negative value is possible due to a bug: https://github.com/broadinstitute/gatk/issues/5492
        if (readPositionByAllele.get(0) > -1 && readPositionByAllele.get(0) < MTFAC.minMedianReadPosition) {
            filterResult.addFilter(GATKVCFConstants.READ_POSITION_FILTER_NAME);
        }
    }

    private void applyGermlineVariantFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY) && vc.hasAttribute(GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE)) {
            final double[] tumorLog10OddsIfSomatic = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
            final Optional<double[]> normalLods = vc.hasAttribute(GATKVCFConstants.NORMAL_LOD_KEY) ?
                    Optional.of(GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.NORMAL_LOD_KEY)) : Optional.empty();
            final double[] negativeLog10AlleleFrequencies = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE);
            final double[] populationAlleleFrequencies = MathUtils.applyToArray(negativeLog10AlleleFrequencies, x -> Math.pow(10,-x));

            final MutableDouble weightedSumOfMafs = new MutableDouble(0);
            for (final Genotype tumorGenotype : vc.getGenotypes()) {
                final String sample = tumorGenotype.getSampleName();
                if (normalSamples.contains(sample)) {
                    continue;
                }

                final List<MinorAlleleFractionRecord> segments = tumorSegments.containsKey(sample) ? tumorSegments.get(sample).getOverlaps(vc).stream().collect(Collectors.toList())
                        : Collections.emptyList();

                // minor allele fraction -- we abbreviate the name to make the formulas below less cumbersome
                final double maf = segments.isEmpty() ? 0.5 : segments.get(0).getMinorAlleleFraction();

                weightedSumOfMafs.add(maf * MathUtils.sum(tumorGenotype.getAD()));
            }

            // note that this includes the ref
            final int[] alleleCounts = sumADsOverSamples(vc, true, false);

            // weighted average of sample minor allele fractions.  This is the expected allele fraction of a germline het in the aggregated read counts
            final double maf = weightedSumOfMafs.getValue() / MathUtils.sum(alleleCounts);

            // exclude the ref
            final int[] altCounts = Arrays.copyOfRange(alleleCounts, 1, alleleCounts.length);

            final int refCount = alleleCounts[0];

            final int totalCount = (int) MathUtils.sum(alleleCounts);

            final double[] altAlleleFractions = Arrays.stream(altCounts)
                    .mapToDouble(c -> c == 0 ? 0 : ((double) c) / totalCount).toArray();

            // this is \chi in the docs, the correction factor for tumor likelihoods if forced to have maf or 1 - maf
            // as the allele fraction
            final double[] log10OddsOfGermlineHetVsSomatic = new IndexRange(0, altAlleleFractions.length).mapToDouble(n -> {
                if (altCounts[n] + refCount == 0) {
                    return 0;
                }
                final double log10GermlineAltMinorLikelihood = log10PowAB(1-maf, refCount) + log10PowAB(maf, altCounts[n]);
                final double log10GermlineAltMajorLikelihood = log10PowAB(maf, refCount) + log10PowAB(1-maf, altCounts[n]);
                final double log10GermlineLikelihood = MathUtils.LOG10_ONE_HALF + MathUtils.log10SumLog10(log10GermlineAltMinorLikelihood, log10GermlineAltMajorLikelihood);

                final double log10SomaticLikelihood = log10PowAB(1 - altAlleleFractions[n], refCount) + log10PowAB(altAlleleFractions[n], altCounts[n]);
                return log10GermlineLikelihood - log10SomaticLikelihood;
            });

            // see docs -- basically the tumor likelihood for a germline hom alt is approximately equal to the somatic likelihood
            // as long as the allele fraction is high
            final double[] log10OddsOfGermlineHomAltVsSomatic = MathUtils.applyToArray(altAlleleFractions, x-> x < MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT ? Double.NEGATIVE_INFINITY : 0);

            final double[] log10GermlinePosteriors = GermlineProbabilityCalculator.calculateGermlineProbabilities(
                    populationAlleleFrequencies, log10OddsOfGermlineHetVsSomatic, log10OddsOfGermlineHomAltVsSomatic, normalLods, MTFAC.log10PriorProbOfSomaticEvent);

            filterResult.addAttribute(GATKVCFConstants.GERMLINE_QUAL_VCF_ATTRIBUTE, Arrays.stream(log10GermlinePosteriors).mapToInt(p -> (int) QualityUtils.phredScaleLog10ErrorRate(p)).toArray());
            final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLog10OddsIfSomatic);
            if (log10GermlinePosteriors[indexOfMaxTumorLod] > Math.log10(MTFAC.maxGermlinePosterior)) {
                filterResult.addFilter(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME);
            }
        }
    }

    private void applyInsufficientEvidenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY)) {
            final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);

            if (MathUtils.arrayMax(tumorLods) < MTFAC.TUMOR_LOD_THRESHOLD) {
                filterResult.addFilter(GATKVCFConstants.TUMOR_LOD_FILTER_NAME);
            }
        }
    }

    // filter out anything called in tumor that would also be called in the normal if it were treated as a tumor.
    // this handles shared artifacts, such as ones due to alignment and any shared aspects of sequencing
    private void applyArtifactInNormalFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (!( vc.hasAttribute(GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE)
                && vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY))) {
            return;
        }

        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        final int[] tumorAlleleDepths = sumADsOverSamples(vc, true, false);
        final int tumorDepth = (int) MathUtils.sum(tumorAlleleDepths);
        final int tumorAltDepth = tumorAlleleDepths[indexOfMaxTumorLod + 1];

        final int[] normalAlleleDepths = sumADsOverSamples(vc, false, true);
        final int normalDepth = (int) MathUtils.sum(normalAlleleDepths);
        final int normalAltDepth = normalAlleleDepths[indexOfMaxTumorLod + 1];

        // if normal AF << tumor AF, don't filter regardless of LOD
        final double tumorAlleleFraction = (double) tumorAltDepth / tumorDepth;
        final double normalAlleleFraction = normalDepth == 0 ? 0 : (double) normalAltDepth / normalDepth;

        if (normalAlleleFraction < MIN_NORMAL_ARTIFACT_RATIO * tumorAlleleFraction)  {
            return;
        }

        // negative because vcf shows log odds of not artifact / artifact (in order to have bigger positive --> good)
        final double[] normalArtifactLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE);
        if (-normalArtifactLods[indexOfMaxTumorLod] > MTFAC.NORMAL_ARTIFACT_LOD_THRESHOLD) {
            filterResult.addFilter(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME);
            return;
        }

        // the above filter misses artifacts whose support in the normal consists entirely of low base quality reads
        // Since a lot of low-BQ reads is itself evidence of an artifact, we filter these by hand via an estimated LOD
        // that uses the average base quality of *ref* reads in the normal
        final int medianRefBaseQuality = vc.getAttributeAsIntList(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY, IMPUTED_NORMAL_BASE_QUALITY).get(0);
        final double normalPValue = 1 - new BinomialDistribution(null, normalDepth, QualityUtils.qualToErrorProb(medianRefBaseQuality))
                .cumulativeProbability(normalAltDepth - 1);

        if (normalPValue < M2FiltersArgumentCollection.normalPileupPValueThreshold) {
            filterResult.addFilter(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME);
        }
    }

    private void applyStrandArtifactFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (!vc.hasAttribute(GATKVCFConstants.STRAND_ARTIFACT_POSTERIOR_KEY)) {
            return;
        }

        final List<Double> posteriorProbabilities = vc.getAttributeAsDoubleList(GATKVCFConstants.STRAND_ARTIFACT_POSTERIOR_KEY, 0.0);
        final List<Double> mapAlleleFractionEstimates = vc.getAttributeAsDoubleList(GATKVCFConstants.STRAND_ARTIFACT_AF_KEY, 0.0);

        if (posteriorProbabilities == null || mapAlleleFractionEstimates == null){
            return;
        }

        final int maxZIndex = MathUtils.maxElementIndex(Doubles.toArray(posteriorProbabilities));

        if (maxZIndex == StrandArtifact.ArtifactState.NO_ARTIFACT.ordinal()){
            return;
        }

        if (posteriorProbabilities.get(maxZIndex) > MTFAC.strandArtifactPosteriorProbThreshold &&
                mapAlleleFractionEstimates.get(maxZIndex) < MTFAC.strandArtifactAlleleFractionThreshold){
            filterResult.addFilter(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME);
        }
    }

    private void applyClusteredEventFilter(final VariantContext vc, final FilterResult filterResult) {
        final Integer eventCount = vc.getAttributeAsInt(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, -1);
        if (eventCount > MTFAC.maxEventsInRegion) {
            filterResult.addFilter(GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME);
        }
    }

    // This filter checks for the case in which PCR-duplicates with unique UMIs (which we assume is caused by false adapter priming)
    // amplify the erroneous signal for an alternate allele.
    private void applyDuplicatedAltReadFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if (!vc.hasAttribute(UniqueAltReadCount.KEY)) {
            return;
        }

        final int uniqueReadSetCount = vc.getAttributeAsInt(UniqueAltReadCount.KEY, 1);

        if (uniqueReadSetCount <= MTFAC.uniqueAltReadCount) {
            filterResult.addFilter(GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME);
        }
    }

     private void applyReadOrientationFilter(final VariantContext vc, final FilterResult filterResult, final Optional<FilteringFirstPass> firstPass){
        if (! vc.isSNP()){
            return;
        }

         final List<ImmutablePair<Integer, Double>> depthsAndPosteriors = new ArrayList<>();

         for (final Genotype tumorGenotype : vc.getGenotypes()) {
             if (normalSamples.contains(tumorGenotype.getSampleName())) {
                 continue;
             }

             if (! tumorGenotype.hasExtendedAttribute(GATKVCFConstants.ROF_POSTERIOR_KEY) || ! tumorGenotype.hasExtendedAttribute(GATKVCFConstants.ROF_PRIOR_KEY)){
                 continue;
             }

             final double artifactPosterior = GATKProtectedVariantContextUtils.getAttributeAsDouble(tumorGenotype, GATKVCFConstants.ROF_POSTERIOR_KEY, 0.0);
             final int[] ADs = tumorGenotype.getAD();
             final int altCount = (int) MathUtils.sum(ADs) - ADs[0];

             depthsAndPosteriors.add(ImmutablePair.of(altCount, artifactPosterior));
         }

        final double artifactPosterior = weightedMedianPosteriorProbability(depthsAndPosteriors);

        if (! firstPass.isPresent()) {
            // During first pass we simply collect the posterior artifact probabilities
            filterResult.setReadOrientationPosterior(artifactPosterior);
            return;
        } else {
            final double threshold = firstPass.get().getFilterStats(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME).getThreshold();

            if (artifactPosterior > threshold){
                filterResult.addFilter(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME);
            }
        }
    }

    private void applyFilteredHaplotypeFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult, final Optional<FilteringFirstPass> firstPass){
        if ( firstPass.isPresent() && firstPass.get().isOnFilteredHaplotype(vc, MTFAC.maxDistanceToFilteredCallOnSameHaplotype)){
            filterResult.addFilter(GATKVCFConstants.BAD_HAPLOTYPE_FILTER_NAME);
        }
    }

    private void applyStrictStrandFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {

        if (! MTFAC.strictStrandBias) {
            return;
        }

        final MutableInt altForwardCount = new MutableInt(0);
        final MutableInt altReverseCount = new MutableInt(0);

        for (final Genotype g : vc.getGenotypes()) {
            if (normalSamples.contains(g.getSampleName())) {
                continue;
            } else if (!g.hasExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)) {
                return;
            } else {
                final int[] strandBiasCounts = GATKProtectedVariantContextUtils.getAttributeAsIntArray(g, GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, () -> null, 0);
                altForwardCount.add(StrandBiasBySample.getAltForwardCountFromFlattenedContingencyTable(strandBiasCounts));
                altReverseCount.add(StrandBiasBySample.getAltReverseCountFromFlattenedContingencyTable(strandBiasCounts));
            }
        }

        // filter if there is no alt evidence in the forward or reverse strand
        if ( altForwardCount.getValue() == 0 || altReverseCount.getValue() == 0) {
            filterResult.addFilter(GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME);
        }
    }

    private void applyNRatioFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        final int[] ADs = sumADsOverSamples(vc, true, true);
        final int altCount = (int) MathUtils.sum(ADs) - ADs[0];
      
        // if there is no NCount annotation or the altCount is 0, don't apply the filter
        if (altCount == 0 ) {
            return;
        }

        final int NCount = vc.getAttributeAsInt(GATKVCFConstants.N_COUNT_KEY,0);

        if ((double) NCount / altCount >= MTFAC.nRatio ) {
            filterResult.addFilter(GATKVCFConstants.N_RATIO_FILTER_NAME);
        }
    }

    private void applyChimericOriginalAlignmentFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {

        if (vc.hasAttribute(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY) && vc.isBiallelic()) {
            final int altCount = vc.getGenotypes().stream().mapToInt(g -> g.getAD()[1]).sum();
            final int nonMtOa = vc.getAttributeAsInt(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY, 0);
            if ((double) nonMtOa / altCount > MTFAC.nonMtAltByAlt) {
                filterResult.addFilter(GATKVCFConstants.CHIMERIC_ORIGINAL_ALIGNMENT_FILTER_NAME);
            }
        }
    }
    private void applyLODFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if(vc.isBiallelic()) {
            final Double lod = vc.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, 1);
            final Double depth = vc.getAttributeAsDouble(VCFConstants.DEPTH_KEY, 1);
            final Double lodByDepth = lod / depth;
            if (lodByDepth < MTFAC.lodByDepth) {
                filterResult.addFilter(GATKVCFConstants.LOW_AVG_ALT_QUALITY_FILTER_NAME);
            }
        }
    }

    public FilterResult calculateFilters(final M2FiltersArgumentCollection MTFAC, final VariantContext vc,
                                         final Optional<FilteringFirstPass> firstPass) {
        firstPass.ifPresent(ffp -> Utils.validate(ffp.isReadyForSecondPass(), "First pass information has not been processed into a model for the second pass."));
        final FilterResult filterResult = new FilterResult();
        applyInsufficientEvidenceFilter(MTFAC, vc, filterResult);
        applyDuplicatedAltReadFilter(MTFAC, vc, filterResult);
        applyStrandArtifactFilter(MTFAC, vc, filterResult);
        applyBaseQualityFilter(MTFAC, vc, filterResult);
        applyMappingQualityFilter(MTFAC, vc, filterResult);
        applyContaminationFilter(MTFAC, vc, filterResult);

        if (!MTFAC.mitochondria) {
            applyFilteredHaplotypeFilter(MTFAC, vc, filterResult, firstPass);
            applyClusteredEventFilter(vc, filterResult);
            applyTriallelicFilter(vc, filterResult);
            applyPanelOfNormalsFilter(MTFAC, vc, filterResult);
            applyGermlineVariantFilter(MTFAC, vc, filterResult);
            applyArtifactInNormalFilter(MTFAC, vc, filterResult);
            applySTRFilter(vc, filterResult);
            applyMedianFragmentLengthDifferenceFilter(MTFAC, vc, filterResult);
            applyReadPositionFilter(MTFAC, vc, filterResult);
            // The ReadOrientation filter uses the information gathered during the first pass
            applyReadOrientationFilter(vc, filterResult, firstPass);
            applyStrictStrandFilter(MTFAC, vc, filterResult);
            applyNRatioFilter(MTFAC, vc, filterResult);
        } else {
            applyChimericOriginalAlignmentFilter(MTFAC, vc, filterResult);
            applyLODFilter(MTFAC, vc, filterResult);
        }

        return filterResult;
    }

    // log10(a^b) = b * log10(a) AND if b = a = 0 the result should be 0, not NaN.  This applies when a is a binomial
    // probability of success and b is the success count -- the likelihood is zero
    private static double log10PowAB(final double a, final int b) {
        return b == 0 ? 0 : b * Math.log10(a);
    }
}
