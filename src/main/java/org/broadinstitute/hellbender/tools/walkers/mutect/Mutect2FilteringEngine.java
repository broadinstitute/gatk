package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.State;

import org.broadinstitute.hellbender.tools.walkers.readorientation.Hyperparameters;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by David Benjamin on 9/15/16.
 */
public class Mutect2FilteringEngine {
    public static final double MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT = 0.9;
    private M2FiltersArgumentCollection MTFAC;
    private final double contamination;
    private final double somaticPriorProb;
    private final String tumorSample;
    final OverlapDetector<MinorAlleleFractionRecord> tumorSegments;
    public static final String FILTERING_STATUS_VCF_KEY = "filtering_status";
    private final List<Hyperparameters> hyperparametersForReadOrientaitonModel;

    public Mutect2FilteringEngine(final M2FiltersArgumentCollection MTFAC, final String tumorSample) {
        this.MTFAC = MTFAC;
        contamination = MTFAC.contaminationTable == null ? 0.0 : ContaminationRecord.readFromFile(MTFAC.contaminationTable).get(0).getContamination();
        this.tumorSample = tumorSample;
        somaticPriorProb = Math.pow(10, MTFAC.log10PriorProbOfSomaticEvent);

        final List<MinorAlleleFractionRecord> tumorMinorAlleleFractionRecords = MTFAC.tumorSegmentationTable == null ?
                Collections.emptyList() : MinorAlleleFractionRecord.readFromFile(MTFAC.tumorSegmentationTable);
        tumorSegments = OverlapDetector.create(tumorMinorAlleleFractionRecords);
        hyperparametersForReadOrientaitonModel = Hyperparameters.readHyperparameters(MTFAC.hyperparameterTable);
    }

    private void applyContaminationFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        final Genotype tumorGenotype = vc.getGenotype(tumorSample);
        final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(tumorGenotype, VCFConstants.ALLELE_FREQUENCY_KEY,
                () -> new double[] {1.0}, 1.0);
        final int maxFractionIndex = MathUtils.maxElementIndex(alleleFractions);
        final int[] ADs = tumorGenotype.getAD();
        final int altCount = ADs[maxFractionIndex + 1];   // AD is all alleles, while AF is alts only, hence the +1 offset
        final int depth = (int) MathUtils.sum(ADs);
        final double alleleFrequency = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc,
                GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE, () -> new double[] {0.0, 0.0}, 0)[maxFractionIndex];

        final double somaticLikelihood = 1.0 / (depth + 1);

        final double singleContaminantLikelihood = 2 * alleleFrequency * (1 - alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination/2)
                + MathUtils.square(alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination);
        final double manyContaminantLikelihood = MathUtils.binomialProbability(depth, altCount, contamination * alleleFrequency);
        final double contaminantLikelihood = Math.max(singleContaminantLikelihood, manyContaminantLikelihood);
        final double posteriorProbOfContamination = (1 - somaticPriorProb) * contaminantLikelihood / ((1 - somaticPriorProb) * contaminantLikelihood + somaticPriorProb * somaticLikelihood);

        vcb.attribute(GATKVCFConstants.POSTERIOR_PROB_OF_CONTAMINATION_ATTRIBUTE, posteriorProbOfContamination);
        if (posteriorProbOfContamination > MTFAC.maxContaminationProbability) {
            vcb.filter(GATKVCFConstants.CONTAMINATION_FILTER_NAME);
        }
    }

    private void applyTriallelicFilter(final VariantContext vc, final VariantContextBuilder vcb) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY)) {
            final double[] tumorLods = getDoubleArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
            final long numPassingAltAlleles = Arrays.stream(tumorLods).filter(x -> x > MTFAC.TUMOR_LOD_THRESHOLD).count();

            if (numPassingAltAlleles > MTFAC.numAltAllelesThreshold) {
                vcb.filter(GATKVCFConstants.MULTIALLELIC_FILTER_NAME);
            }
        }
    }

    private static void applySTRFilter(final VariantContext vc, final VariantContextBuilder vcb) {
        // STR contractions, such as ACTACTACT -> ACTACT, are overwhelmingly false positives so we hard filter by default
        if (vc.isIndel()) {
            final int[] rpa = vc.getAttributeAsList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY).stream()
                    .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
            final String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");
            if (rpa != null && rpa.length > 1 && ru.length() > 1) {
                final int refCount = rpa[0];
                final int altCount = rpa[1];

                if (refCount - altCount == 1) {
                    vcb.filter(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME);
                }
            }
        }
    }

    private static void applyPanelOfNormalsFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        final boolean siteInPoN = vc.hasAttribute(GATKVCFConstants.IN_PON_VCF_ATTRIBUTE);
        if (siteInPoN) {
            vcb.filter(GATKVCFConstants.PON_FILTER_NAME);
        }
    }

    private void applyMedianBaseQualityDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        final int[] baseQualityByAllele = getIntArrayTumorField(vc, BaseQuality.KEY);
        if (baseQualityByAllele != null && baseQualityByAllele[0] < MTFAC.minMedianBaseQuality) {
            vcb.filter(GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME);
        }
    }

    private void applyMedianMappingQualityDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        final int[] mappingQualityByAllele = getIntArrayTumorField(vc, MappingQuality.KEY);
        if (mappingQualityByAllele != null && mappingQualityByAllele[0] < MTFAC.minMedianMappingQuality) {
            vcb.filter(GATKVCFConstants.MEDIAN_MAPPING_QUALITY_FILTER_NAME);
        }
    }

    private void applyMedianFragmentLengthDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        final int[] fragmentLengthByAllele = getIntArrayTumorField(vc, FragmentLength.KEY);
        if (fragmentLengthByAllele != null && Math.abs(fragmentLengthByAllele[1] - fragmentLengthByAllele[0]) > MTFAC.maxMedianFragmentLengthDifference) {
            vcb.filter(GATKVCFConstants.MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME);
        }
    }

    private void applyReadPositionFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        final int[] readPositionByAllele = getIntArrayTumorField(vc, ReadPosition.KEY);
        if (readPositionByAllele != null) {
            final int insertionSize =  Math.max(vc.getAltAlleleWithHighestAlleleCount().getBases().length - vc.getReference().getBases().length, 0);
            if (insertionSize + readPositionByAllele[0] < MTFAC.minMedianReadPosition) {
                vcb.filter(GATKVCFConstants.READ_POSITION_FILTER_NAME);
            }
        }
    }

    private void applyGermlineVariantFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY) && vc.hasAttribute(GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE)) {
            final double[] tumorLog10OddsIfSomatic = getDoubleArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
            final Optional<double[]> normalLods = vc.hasAttribute(GATKVCFConstants.NORMAL_LOD_KEY) ?
                    Optional.of(getDoubleArrayAttribute(vc, GATKVCFConstants.NORMAL_LOD_KEY)) : Optional.empty();
            final double[] populationAlleleFrequencies = getDoubleArrayAttribute(vc, GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE);

            final List<MinorAlleleFractionRecord> segments = tumorSegments.getOverlaps(vc).stream().collect(Collectors.toList());

            // minor allele fraction -- we abbreviate the name to make the formulas below less cumbersome
            final double maf = segments.isEmpty() ? 0.5 : segments.get(0).getMinorAlleleFraction();

            final double[] altAlleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc.getGenotype(tumorSample), GATKVCFConstants.ALLELE_FRACTION_KEY, () -> null, 0);

            // note that this includes the ref
            final int[] alleleCounts = vc.getGenotype(tumorSample).getAD();
            // exclude the ref
            final int[] altCounts = Arrays.copyOfRange(alleleCounts, 1, alleleCounts.length);

            final int refCount = alleleCounts[0];

            // this is \chi in the docs, the correction factor for tumor likelihoods if forced to have maf or 1 - maf
            // as the allele fraction
            final double[] log10OddsOfGermlineHetVsSomatic = new IndexRange(0, altAlleleFractions.length).mapToDouble(n -> {
                final double log10GermlineAltMinorLikelihood = refCount * Math.log10(1 - maf) + altCounts[n] * Math.log10(maf);
                final double log10GermlineAltMajorLikelihood = refCount * Math.log10(maf) + altCounts[n] * Math.log10(1 - maf);
                final double log10GermlineLikelihood = MathUtils.LOG10_ONE_HALF + MathUtils.log10SumLog10(log10GermlineAltMinorLikelihood, log10GermlineAltMajorLikelihood);

                final double log10SomaticLikelihood = refCount * Math.log10(1 - altAlleleFractions[n]) + altCounts[n] * Math.log10(altAlleleFractions[n]);
                return log10GermlineLikelihood - log10SomaticLikelihood;
            });

            // see docs -- basically the tumor likelihood for a germline hom alt is approximately equal to the somatic likelihood
            // as long as the allele fraction is high
            final double[] log10OddsOfGermlineHomAltVsSomatic = MathUtils.applyToArray(altAlleleFractions, x-> x < MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT ? Double.NEGATIVE_INFINITY : 0);

            final double[] log10GermlinePosteriors = GermlineProbabilityCalculator.calculateGermlineProbabilities(
                    populationAlleleFrequencies, log10OddsOfGermlineHetVsSomatic, log10OddsOfGermlineHomAltVsSomatic, normalLods, MTFAC.log10PriorProbOfSomaticEvent);

            vcb.attribute(GATKVCFConstants.GERMLINE_POSTERIORS_VCF_ATTRIBUTE, log10GermlinePosteriors);
            final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLog10OddsIfSomatic);
            if (log10GermlinePosteriors[indexOfMaxTumorLod] > Math.log10(MTFAC.maxGermlinePosterior)) {
                vcb.filter(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME);
            }
        }
    }

    private static void applyInsufficientEvidenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY)) {
            final double[] tumorLods = getDoubleArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);

            if (MathUtils.arrayMax(tumorLods) < MTFAC.TUMOR_LOD_THRESHOLD) {
                vcb.filter(GATKVCFConstants.TUMOR_LOD_FILTER_NAME);
            }
        }
    }

    // filter out anything called in tumor that would also be called in the normal if it were treated as a tumor.
    // this handles shared artifacts, such as ones due to alignment and any shared aspects of sequencing
    private static void applyArtifactInNormalFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        if (!( vc.hasAttribute(GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE)
                && vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY))) {
            return;
        }

        final double[] normalArtifactLods = getDoubleArrayAttribute(vc, GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE);
        final double[] tumorLods = getDoubleArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        if (normalArtifactLods[indexOfMaxTumorLod] > MTFAC.NORMAL_ARTIFACT_LOD_THRESHOLD) {
            vcb.filter(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME);
        }
    }

    private static double[] getDoubleArrayAttribute(final VariantContext vc, final String attribute) {
        return GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, attribute, () -> null, -1);
    }

    private void applyStrandArtifactFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        Genotype tumorGenotype = vc.getGenotype(tumorSample);
        final double[] posteriorProbabilities = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(
                tumorGenotype, (GATKVCFConstants.POSTERIOR_PROBABILITIES_KEY), () -> null, -1);
        final double[] mapAlleleFractionEstimates = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(
                tumorGenotype, (GATKVCFConstants.MAP_ALLELE_FRACTIONS_KEY), () -> null, -1);

        if (posteriorProbabilities == null || mapAlleleFractionEstimates == null){
            return;
        }

        final int maxZIndex = MathUtils.maxElementIndex(posteriorProbabilities);

        if (maxZIndex == StrandArtifact.ArtifactState.NO_ARTIFACT.ordinal()){
            return;
        }

        if (posteriorProbabilities[maxZIndex] > MTFAC.strandArtifactPosteriorProbThreshold &&
                mapAlleleFractionEstimates[maxZIndex] < MTFAC.strandArtifactAlleleFractionThreshold){
            vcb.filter(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME);
        }
    }

    private void applyClusteredEventFilter(final VariantContext vc, final VariantContextBuilder vcb) {
        final Integer eventCount = vc.getAttributeAsInt(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, -1);
        if (eventCount > MTFAC.maxEventsInRegion) {
            vcb.filter(GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME);
        }
    }

    // This filter checks for the case in which PCR-duplicates with unique UMIs (which we assume is caused by false adapter priming)
    // amplify the erroneous signal for an alternate allele.
    private void applyDuplicatedAltReadFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        Genotype tumorGenotype = vc.getGenotype(tumorSample);

        if (!tumorGenotype.hasExtendedAttribute(UniqueAltReadCount.UNIQUE_ALT_READ_SET_COUNT_KEY)) {
            return;
        }

        final int uniqueReadSetCount = GATKProtectedVariantContextUtils.getAttributeAsInt(tumorGenotype, UniqueAltReadCount.UNIQUE_ALT_READ_SET_COUNT_KEY, -1);

        if (uniqueReadSetCount <= MTFAC.uniqueAltReadCount) {
            vcb.filter(GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME);
        }
    }


    private void applyReadOrientationFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters){
        /***
         * This filter needs {@code REFERENCE_CONTEXT_KEY} and {@code F1R2_KEY}
         */

        // skip INDELs
        if (vc.getAlternateAllele(0).length() > 1){
            return;
        }

        final String referenceContext = vc.getAttributeAsString(GATKVCFConstants.REFERENCE_CONTEXT_KEY, "");
        if (referenceContext.isEmpty()){
            return;
        }

        Utils.validate(referenceContext.length() == 3, String.format("reference context must be of length 3 but got %s", referenceContext));

        // VCF must contain: reference context
        final Nucleotide altAllele = Nucleotide.valueOf(vc.getAlternateAllele(0).toString());
        final State[] relevantArtifactStates;
        switch (altAllele) {
            case A : relevantArtifactStates = new State[]{ State.F1R2_A, State.F2R1_A }; break;
            case C : relevantArtifactStates = new State[]{ State.F1R2_C, State.F2R1_C }; break;
            case G : relevantArtifactStates = new State[]{ State.F1R2_G, State.F2R1_G }; break;
            case T : relevantArtifactStates = new State[]{ State.F1R2_T, State.F2R1_T }; break;
            default: throw new UserException(String.format("Alt allele must be in {A, C, G, T} but got %s", altAllele));
        }

        final Genotype tumorGenotype = vc.getGenotype(tumorSample);
        final int depth = tumorGenotype.getDP();
        final int altDepth = tumorGenotype.getAD()[1]; // what about multiple alleles?
        final int altF1R2Depth = vc.getAttributeAsIntList(GATKVCFConstants.F1R2_KEY, -1).get(1); // f1r2 genotype field, vc.getAttributeAsIntList()
        final Hyperparameters hyps = hyperparametersForReadOrientaitonModel.stream()
                .filter(h -> h.getReferenceContext().equals(referenceContext))
                .findFirst()
                .get();

        // A by K matrix of prior probabilities over K latent states, given allele a \in A
        // \pi_{ak} is the prior probability of state k given observed allele a.
        final double[] pi = hyps.getPi();

        // a vector of length K, the probability of drawing an alt read (i.e. allele fraction) given z
        final double[] f = hyps.getF();

        // a vector of length K, the probability of drawing an F1R2 alt read given z
        final double[] theta = hyps.getTheta();

        final double[] log10UnnormalizedPosteriorProbabilities = new double[LearnHyperparametersEngine.NUM_STATES];
        for (State z : State.values()){
            final int k = z.ordinal();
            if (z == State.SOMATIC_HET){
                log10UnnormalizedPosteriorProbabilities[k] = Math.log10(pi[k]) +
                        MathUtils.log10BetaBinomialDensity(altDepth, depth, LearnHyperparametersEngine.alpha, LearnHyperparametersEngine.beta) +
                        MathUtils.log10BinomialProbability(altDepth, altF1R2Depth, Math.log10(theta[k]));

            } else if (State.getNonArtifactStates().contains(z)) {
                log10UnnormalizedPosteriorProbabilities[k] = Math.log10(pi[k]) +
                        MathUtils.log10BinomialProbability(depth, altDepth, f[k]) +
                        MathUtils.log10BinomialProbability(altDepth, altF1R2Depth, theta[k]);

            } else {
                // we are in artifact states
                if (altAllele != State.getAltAlleleOfTransition(z)){
                    // We assume that the artifact states that are irrelevant to the particular transition we're looking at
                    // did not happen e.g. under G -> T SNP, F1R2_A, F1R2_C, and F1R2_G are all irrelevant
                    // In other words, we compute the posterior probability of F1R2_T and F2R1_T given that these irrelevant
                    // states did not occur. To do so we set the log posterior probabilities of these states to -Infinity
                    log10UnnormalizedPosteriorProbabilities[k] = Double.NEGATIVE_INFINITY;
                } else {
                    log10UnnormalizedPosteriorProbabilities[k] = Math.log10(pi[k]) +
                            MathUtils.log10BinomialProbability(depth, altDepth, f[k]) +
                            MathUtils.log10BinomialProbability(altDepth, altF1R2Depth, theta[k]);
                }

            }

        }

        final double[] posteriorProbabilties = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedPosteriorProbabilities);

        for (State z : relevantArtifactStates){
            if (posteriorProbabilties[z.ordinal()] > MTFAC.readOrientationFilterThreshold){
                filters.add(GATKVCFConstants.READ_ORIENTATION_FILTER_NAME);
            }
        }
    }


    //TODO: building a list via repeated side effects is ugly
    public void applyFilters(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final VariantContextBuilder vcb) {
        vcb.filters(new HashSet<>());
        applyInsufficientEvidenceFilter(MTFAC, vc, vcb);
        applyClusteredEventFilter(vc, vcb);
        applyDuplicatedAltReadFilter(MTFAC, vc, vcb);
        applyTriallelicFilter(vc, vcb);
        applyPanelOfNormalsFilter(MTFAC, vc, vcb);
        applyGermlineVariantFilter(MTFAC, vc, vcb);
        applyArtifactInNormalFilter(MTFAC, vc, vcb);
        applyStrandArtifactFilter(MTFAC, vc, vcb);
        applySTRFilter(vc, vcb);
        applyContaminationFilter(MTFAC, vc, vcb);
        applyMedianBaseQualityDifferenceFilter(MTFAC, vc, vcb);
        applyMedianMappingQualityDifferenceFilter(MTFAC, vc, vcb);
        applyMedianFragmentLengthDifferenceFilter(MTFAC, vc, vcb);
        applyReadPositionFilter(MTFAC, vc, vcb);
        applyReadOrientationFilter(MTFAC, vc, vcb);
    }

    private int[] getIntArrayTumorField(final VariantContext vc, final String key) {
        return GATKProtectedVariantContextUtils.getAttributeAsIntArray(vc.getGenotype(tumorSample), key, () -> null, 0);
    }
}
