package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.tools.walkers.readorientation.Hyperparameters;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine;
import org.broadinstitute.hellbender.tools.walkers.readorientation.LearnHyperparametersEngine.State;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Created by David Benjamin on 9/15/16.
 */
public class Mutect2FilteringEngine {
    private M2FiltersArgumentCollection MTFAC;
    private final double contamination;
    private final String tumorSample;
    public static final String FILTERING_STATUS_VCF_KEY = "filtering_status";

    public Mutect2FilteringEngine(final M2FiltersArgumentCollection MTFAC, final String tumorSample) {
        this.MTFAC = MTFAC;
        contamination = MTFAC.contaminationTable == null ? 0.0 : ContaminationRecord.readContaminationTable(MTFAC.contaminationTable).get(0).getContamination();
        this.tumorSample = tumorSample;
    }

    // very naive M1-style contamination filter -- remove calls with AF less than the contamination fraction
    private void applyContaminationFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final Genotype tumorGenotype = vc.getGenotype(tumorSample);
        final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(tumorGenotype, VCFConstants.ALLELE_FREQUENCY_KEY,
                () -> new double[] {1.0}, 1.0);
        final double maxFraction = MathUtils.arrayMax(alleleFractions);
        if (maxFraction < contamination) {
            filters.add(GATKVCFConstants.CONTAMINATION_FILTER_NAME);
        }
    }

    private void applyTriallelicFilter(final VariantContext vc, final Collection<String> filters) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY)) {
            final double[] tumorLods = getArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
            final long numPassingAltAlleles = Arrays.stream(tumorLods).filter(x -> x > MTFAC.TUMOR_LOD_THRESHOLD).count();

            if (numPassingAltAlleles > MTFAC.numAltAllelesThreshold) {
                filters.add(GATKVCFConstants.MULTIALLELIC_FILTER_NAME);
            }
        }
    }

    private static void applySTRFilter(final VariantContext vc, final Collection<String> filters) {
        // STR contractions, such as ACTACTACT -> ACTACT, are overwhelmingly false positives so we hard filter by default
        if (vc.isIndel()) {
            final int[] rpa = vc.getAttributeAsList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY).stream()
                    .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
            final String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");
            if (rpa != null && rpa.length > 1 && ru.length() > 1) {
                final int refCount = rpa[0];
                final int altCount = rpa[1];

                if (refCount - altCount == 1) {
                    filters.add(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME);
                }
            }
        }
    }

    private static void applyPanelOfNormalsFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final boolean siteInPoN = vc.hasAttribute(GATKVCFConstants.IN_PON_VCF_ATTRIBUTE);
        if (siteInPoN) {
            filters.add(GATKVCFConstants.PON_FILTER_NAME);
        }
    }

    private void applyMedianBaseQualityDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] baseQualityByAllele = getIntArrayTumorField(vc, BaseQuality.KEY);
        if (baseQualityByAllele != null && baseQualityByAllele[0] < MTFAC.minMedianBaseQuality) {
            filters.add(GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME);
        }
    }

    private void applyMedianMappingQualityDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] mappingQualityByAllele = getIntArrayTumorField(vc, MappingQuality.KEY);
        if (mappingQualityByAllele != null && mappingQualityByAllele[0] < MTFAC.minMedianMappingQuality) {
            filters.add(GATKVCFConstants.MEDIAN_MAPPING_QUALITY_FILTER_NAME);
        }
    }

    private void applyMedianFragmentLengthDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] fragmentLengthByAllele = getIntArrayTumorField(vc, FragmentLength.KEY);
        if (fragmentLengthByAllele != null && Math.abs(fragmentLengthByAllele[1] - fragmentLengthByAllele[0]) > MTFAC.maxMedianFragmentLengthDifference) {
            filters.add(GATKVCFConstants.MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME);
        }
    }

    private void applyReadPositionFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] readPositionByAllele = getIntArrayTumorField(vc, ReadPosition.KEY);
        if (readPositionByAllele != null) {
            final int insertionSize =  Math.max(vc.getAltAlleleWithHighestAlleleCount().getBases().length - vc.getReference().getBases().length, 0);
            if (insertionSize + readPositionByAllele[0] < MTFAC.minMedianReadPosition) {
                filters.add(GATKVCFConstants.READ_POSITION_FILTER_NAME);
            }
        }
    }



    private static void applyGermlineVariantFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY) && vc.hasAttribute(GATKVCFConstants.GERMLINE_POSTERIORS_VCF_ATTRIBUTE)) {
            final double[] tumorLods = getArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
            final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

            final double[] log10GermlinePosteriors = getArrayAttribute(vc, GATKVCFConstants.GERMLINE_POSTERIORS_VCF_ATTRIBUTE);
            if (log10GermlinePosteriors[indexOfMaxTumorLod] > Math.log10(MTFAC.maxGermlinePosterior)) {
                filters.add(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME);
            }
        }
    }

    private static void applyInsufficientEvidenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        if (vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY)) {
            final double[] tumorLods = getArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);

            if (MathUtils.arrayMax(tumorLods) < MTFAC.TUMOR_LOD_THRESHOLD) {
                filters.add(GATKVCFConstants.TUMOR_LOD_FILTER_NAME);
            }
        }
    }

    // filter out anything called in tumor that would also be called in the normal if it were treated as a tumor.
    // this handles shared artifacts, such as ones due to alignment and any shared aspects of sequencing
    private static void applyArtifactInNormalFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        if (!( vc.hasAttribute(GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE)
                && vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY))) {
            return;
        }

        final double[] normalArtifactLods = getArrayAttribute(vc, GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE);
        final double[] tumorLods = getArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        if (normalArtifactLods[indexOfMaxTumorLod] > MTFAC.NORMAL_ARTIFACT_LOD_THRESHOLD) {
            filters.add(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME);
        }
    }

    private static double[] getArrayAttribute(final VariantContext vc, final String attribute) {
        return GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, attribute, () -> null, -1);
    }

    private void applyStrandArtifactFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        Genotype tumorGenotype = vc.getGenotype(tumorSample);
        final double[] posteriorProbabilities = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(
                tumorGenotype, (GATKVCFConstants.POSTERIOR_PROBABILITIES_KEY), () -> null, -1);
        final double[] mapAlleleFractionEstimates = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(
                tumorGenotype, (GATKVCFConstants.MAP_ALLELE_FRACTIONS_KEY), () -> null, -1);

        if (posteriorProbabilities == null || mapAlleleFractionEstimates == null){
            return;
        }

        final int maxZIndex = MathUtils.maxElementIndex(posteriorProbabilities);

        if (maxZIndex == StrandArtifact.StrandArtifactZ.NO_ARTIFACT.ordinal()){
            return;
        }

        if (posteriorProbabilities[maxZIndex] > MTFAC.strandArtifactPosteriorProbThreshold &&
                mapAlleleFractionEstimates[maxZIndex] < MTFAC.strandArtifactAlleleFractionThreshold){
            filters.add(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME);
        }
    }

    private void applyClusteredEventFilter(final VariantContext vc, final Collection<String> filters) {
        final Integer eventCount = vc.getAttributeAsInt(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, -1);
        if (eventCount > MTFAC.maxEventsInRegion) {
            filters.add(GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME);
        }
    }

    // This filter checks for the case in which PCR-duplicates with unique UMIs (which we assume is caused by false adapter priming)
    // amplify the erroneous signal for an alternate allele.
    private void applyDuplicatedAltReadFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        Genotype tumorGenotype = vc.getGenotype(tumorSample);

        if (!tumorGenotype.hasExtendedAttribute(UniqueAltReadCount.UNIQUE_ALT_READ_SET_COUNT_KEY)) {
            return;
        }

        final int uniqueReadSetCount = GATKProtectedVariantContextUtils.getAttributeAsInt(tumorGenotype, UniqueAltReadCount.UNIQUE_ALT_READ_SET_COUNT_KEY, -1);

        if (uniqueReadSetCount <= MTFAC.uniqueAltReadCount) {
            filters.add(GATKVCFConstants.DUPLICATED_EVIDENCE_FILTER_NAME);
        }
    }

    /***
     * This filter requires the INFO field annotation {@code REFERENCE_CONTEXT_KEY} and {@code F1R2_KEY}
     */
    private void applyReadOrientationFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters){
        final int altAlleleIndex = 1;
        final String referenceContext = vc.getAttributeAsString(GATKVCFConstants.REFERENCE_CONTEXT_KEY, "");
        final int[] f1r2CountsByAllele = GATKProtectedVariantContextUtils.getAttributeAsIntArray(vc.getGenotype(tumorSample),
                GATKVCFConstants.F1R2_KEY, () -> null, -1);
        if (referenceContext.isEmpty() || f1r2CountsByAllele == null){
            return;
        }

        final List<Hyperparameters> hyperparametersForReadOrientaitonModel =
                Hyperparameters.readHyperparameters(MTFAC.hyperparameterTable);

        // skip INDELs
        if (vc.getAlternateAllele(0).length() > 1){
            return;
        }

        Utils.validate(referenceContext.length() == 3, String.format("reference context must be of length 3 but got %s", referenceContext));

        final Nucleotide refAllele = Nucleotide.valueOf(referenceContext.substring(1,1));
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
        final int depth = (int) MathUtils.sum(tumorGenotype.getAD());
        final int altDepth = tumorGenotype.getAD()[altAlleleIndex]; // what about multiple alleles?
        final int f1r2AltCount = f1r2CountsByAllele[altAlleleIndex];
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
                m_nk = 0;
                x_nk = 0;
            }

            log10UnnormalizedPosteriorProbabilities[k] =
                    LearnHyperparametersEngine.computePosterior(z, m_nk, x_nk, r, pi[k], f[k], theta[k]);
        }

        final double[] posteriorProbabilties = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedPosteriorProbabilities);

        for (State z : relevantArtifactStates){
            if (posteriorProbabilties[z.ordinal()] > MTFAC.readOrientationFilterThreshold){
                filters.add(GATKVCFConstants.READ_ORIENTATION_FILTER_NAME);
            }
        }
    }

    //TODO: building a list via repeated side effects is ugly
    public Set<String> calculateFilters(final M2FiltersArgumentCollection MTFAC, final VariantContext vc) {
        final Set<String> filters = new HashSet<>();
        applyInsufficientEvidenceFilter(MTFAC, vc, filters);
        applyClusteredEventFilter(vc, filters);
        applyDuplicatedAltReadFilter(MTFAC, vc, filters);
        applyTriallelicFilter(vc, filters);
        applyPanelOfNormalsFilter(MTFAC, vc, filters);
        applyGermlineVariantFilter(MTFAC, vc, filters);
        applyArtifactInNormalFilter(MTFAC, vc, filters);
        applyStrandArtifactFilter(MTFAC, vc, filters);
        applySTRFilter(vc, filters);
        applyContaminationFilter(MTFAC, vc, filters);
        applyMedianBaseQualityDifferenceFilter(MTFAC, vc, filters);
        applyMedianMappingQualityDifferenceFilter(MTFAC, vc, filters);
        applyMedianFragmentLengthDifferenceFilter(MTFAC, vc, filters);
        applyReadPositionFilter(MTFAC, vc, filters);
        applyReadOrientationFilter(MTFAC, vc, filters);

        return filters;
    }

    private int[] getIntArrayTumorField(final VariantContext vc, final String key) {
        return GATKProtectedVariantContextUtils.getAttributeAsIntArray(vc.getGenotype(tumorSample), key, () -> null, 0);
    }

}
