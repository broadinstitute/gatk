package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.math3.util.DoubleArray;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Created by David Benjamin on 9/15/16.
 */
public class Mutect2FilteringEngine {

    public final static String ARTIFACT_IN_NORMAL_FILTER_NAME = "artifact_in_normal";
    public final static String MEDIAN_BASE_QUALITY_DIFFERENCE_FILTER_NAME = "base_quality";
    public final static String MEDIAN_MAPPING_QUALITY_DIFFERENCE_FILTER_NAME = "mapping_quality";
    public final static String MEDIAN_CLIPPING_DIFFERENCE_FILTER_NAME = "clipping";
    public final static String MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME = "fragment_length";
    public final static String READ_POSITION_FILTER_NAME = "read_position";
    public final static String CONTAMINATION_FILTER_NAME = "contamination";

    public static final List<String> M_2_FILTER_NAMES = Arrays.asList(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME, GATKVCFConstants.PON_FILTER_NAME,
            GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME, GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME,
            GATKVCFConstants.TUMOR_LOD_FILTER_NAME, GATKVCFConstants.GERMLINE_RISK_FILTER_NAME, GATKVCFConstants.TRIALLELIC_SITE_FILTER_NAME,
            GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME);

    private M2FiltersArgumentCollection MTFAC;
    private final double contamination;
    private final String tumorSample;


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
            filters.add(CONTAMINATION_FILTER_NAME);
        }
    }

    private void applyTriallelicFilter(final VariantContext vc, final Collection<String> filters) {
        if (vc.getNAlleles() > (MTFAC.numAltAllelesThreshold + 1)) {
            filters.add(GATKVCFConstants.TRIALLELIC_SITE_FILTER_NAME);
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
        final boolean siteInPoN = vc.hasAttribute(SomaticGenotypingEngine.IN_PON_VCF_ATTRIBUTE);
        if (siteInPoN) {
            filters.add(GATKVCFConstants.PON_FILTER_NAME);
        }
    }

    private void applyMedianBaseQualityDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] baseQualityByAllele = getIntArrayTumorField(vc, BaseQuality.KEY);
        if (baseQualityByAllele != null && baseQualityByAllele[0] - baseQualityByAllele[1] > MTFAC.maxMedianBaseQualityDifference) {
            filters.add(MEDIAN_CLIPPING_DIFFERENCE_FILTER_NAME);
        }
    }

    private void applyMedianMappingQualityDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] mappingQualityByAllele = getIntArrayTumorField(vc, MappingQuality.KEY);
        if (mappingQualityByAllele != null && mappingQualityByAllele[0] - mappingQualityByAllele[1] > MTFAC.maxMedianMappingQualityDifference) {
            filters.add(MEDIAN_CLIPPING_DIFFERENCE_FILTER_NAME);
        }
    }

    private void applyMedianClippingDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] clippingCountsByAllele = getIntArrayTumorField(vc, ClippedBases.KEY);
        if (clippingCountsByAllele != null && clippingCountsByAllele[1] - clippingCountsByAllele[0] > MTFAC.maxMedianClippingDifference) {
            filters.add(MEDIAN_CLIPPING_DIFFERENCE_FILTER_NAME);
        }
    }

    private void applyMedianFragmentLengthDifferenceFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] fragmentLengthByAllele = getIntArrayTumorField(vc, FragmentLength.KEY);
        if (fragmentLengthByAllele != null && Math.abs(fragmentLengthByAllele[1] - fragmentLengthByAllele[0]) > MTFAC.maxMedianFragmentLengthDifference) {
            filters.add(MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME);
        }
    }

    private void applyReadPositionFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        final int[] readPositionByAllele = getIntArrayTumorField(vc, ReadPosition.KEY);
        if (readPositionByAllele != null) {
            final int insertionSize =  Math.max(vc.getAltAlleleWithHighestAlleleCount().getBases().length - vc.getReference().getBases().length, 0);
            if (insertionSize + readPositionByAllele[1] < MTFAC.minMedianReadPosition) {
                filters.add(READ_POSITION_FILTER_NAME);
            }
        }

        // since read position is measured relative to the reference, this can unfairly penalize long indels
        // thus we add the insertion length
        final int insertionSize =  Math.max(vc.getAltAlleleWithHighestAlleleCount().getBases().length - vc.getReference().getBases().length, 0);
        if (insertionSize + vc.getAttributeAsInt(ReadPosition.KEY, 100) < MTFAC.minMedianReadPosition) {
            filters.add(READ_POSITION_FILTER_NAME);
        }
    }



    private static void applyGermlineVariantFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        if (!vc.hasAttribute(GATKVCFConstants.NORMAL_LOD_KEY) || !vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY)) {
            return;
        }

        final double[] tumorLods = getArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        final boolean siteInCosmic = vc.hasAttribute(SomaticGenotypingEngine.IN_COSMIC_VCF_ATTRIBUTE);
        final boolean siteInDbsnp = vc.hasAttribute(SomaticGenotypingEngine.IN_DBSNP_VCF_ATTRIBUTE);
        if (siteInDbsnp && !siteInCosmic ) {
            // take the normal LOD of the best somatic alt allele
            final double normalLod = getArrayAttribute(vc, GATKVCFConstants.NORMAL_LOD_KEY)[indexOfMaxTumorLod];
            if (normalLod < MTFAC.NORMAL_DBSNP_LOD_THRESHOLD) {
                filters.add(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME);
            }
        }
    }

    // filter out anything called in tumor that would also be called in the normal if it were treated as a tumor.
    // this handles shared artifacts, such as ones due to alignment and any shared aspects of sequencing
    private static void applyArtifactInNormalFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        if (!( vc.hasAttribute(SomaticGenotypingEngine.NORMAL_ARTIFACT_LOD_ATTRIBUTE)
                && vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY))) {
            return;
        }

        final double[] normalArtifactLods = getArrayAttribute(vc, SomaticGenotypingEngine.NORMAL_ARTIFACT_LOD_ATTRIBUTE);
        final double[] tumorLods = getArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        if (normalArtifactLods[indexOfMaxTumorLod] > MTFAC.NORMAL_ARTIFACT_LOD_THRESHOLD) {
            filters.add(ARTIFACT_IN_NORMAL_FILTER_NAME);
        }
    }

    private static double[] getArrayAttribute(final VariantContext vc, final String attribute) {
        return GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, attribute, () -> null, -1);
    }

    private static void applyStrandBiasFilter(final M2FiltersArgumentCollection MTFAC, final VariantContext vc, final Collection<String> filters) {
        if (MTFAC.ENABLE_STRAND_ARTIFACT_FILTER) {
            if (vc.hasAttribute(GATKVCFConstants.TLOD_FWD_KEY) && vc.hasAttribute(GATKVCFConstants.TLOD_REV_KEY)
                    && vc.hasAttribute(GATKVCFConstants.TUMOR_SB_POWER_FWD_KEY) && vc.hasAttribute(GATKVCFConstants.TUMOR_SB_POWER_REV_KEY)) {
                final double forwardLod = vc.getAttributeAsDouble(GATKVCFConstants.TLOD_FWD_KEY, 0.0);
                final double reverseLod = vc.getAttributeAsDouble(GATKVCFConstants.TLOD_REV_KEY, 0.0);
                final double forwardPower = vc.getAttributeAsDouble(GATKVCFConstants.TUMOR_SB_POWER_FWD_KEY, 0.0);
                final double reversePower = vc.getAttributeAsDouble(GATKVCFConstants.TUMOR_SB_POWER_REV_KEY, 0.0);
                if ((forwardPower > MTFAC.STRAND_ARTIFACT_POWER_THRESHOLD && forwardLod < MTFAC.STRAND_ARTIFACT_LOD_THRESHOLD) ||
                        (reversePower > MTFAC.STRAND_ARTIFACT_POWER_THRESHOLD && reverseLod < MTFAC.STRAND_ARTIFACT_LOD_THRESHOLD)) {
                    filters.add(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME);
                }
            }
        }
    }

    private void applyEventDistanceFilters(final VariantContext vc, final Collection<String> filters) {
        final Integer eventCount = vc.getAttributeAsInt(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, -1);
        if (eventCount > MTFAC.maxEventsInHaplotype) {
            filters.add(GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME);
        }
    }

    //TODO: building a list via repeated side effects is ugly
    public Set<String> calculateFilters(final M2FiltersArgumentCollection MTFAC, final VariantContext vc) {
        final Set<String> filters = new HashSet<>();
        applyEventDistanceFilters(vc, filters);
        applyTriallelicFilter(vc, filters);
        applyPanelOfNormalsFilter(MTFAC, vc, filters);
        applyGermlineVariantFilter(MTFAC, vc, filters);
        applyArtifactInNormalFilter(MTFAC, vc, filters);
        applyStrandBiasFilter(MTFAC, vc, filters);
        applySTRFilter(vc, filters);
        applyContaminationFilter(MTFAC, vc, filters);
        applyMedianBaseQualityDifferenceFilter(MTFAC, vc, filters);
        applyMedianMappingQualityDifferenceFilter(MTFAC, vc, filters);
        applyMedianClippingDifferenceFilter(MTFAC, vc, filters);
        applyMedianFragmentLengthDifferenceFilter(MTFAC, vc, filters);
        applyReadPositionFilter(MTFAC, vc, filters);

        return filters;
    }

    private int[] getIntArrayTumorField(final VariantContext vc, final String key) {
        return GATKProtectedVariantContextUtils.getAttributeAsIntArray(vc.getGenotype(tumorSample), key, () -> null, 0);
    }

}
