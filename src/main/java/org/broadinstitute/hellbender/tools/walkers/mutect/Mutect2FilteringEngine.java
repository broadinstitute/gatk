package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Created by davidben on 9/15/16.
 */
public class Mutect2FilteringEngine {

    private Mutect2FilteringEngine() { }

    private static void applyTriallelicFilter(final VariantContext vc, final Collection<String> filters) {
        if (vc.getNAlleles() > 2) {
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

    private static void applyClusteredReadPositionFilter(final M2ArgumentCollection MTAC, final VariantContext vc, final Collection<String> filters) {
        if (MTAC.ENABLE_CLUSTERED_READ_POSITION_FILTER) {
            final Double tumorFwdPosMedian = (Double) vc.getAttribute(GATKVCFConstants.MEDIAN_LEFT_OFFSET_KEY);
            final Double tumorRevPosMedian = (Double) vc.getAttribute(GATKVCFConstants.MEDIAN_RIGHT_OFFSET_KEY);
            final Double tumorFwdPosMAD = (Double) vc.getAttribute(GATKVCFConstants.MAD_MEDIAN_LEFT_OFFSET_KEY);
            final Double tumorRevPosMAD = (Double) vc.getAttribute(GATKVCFConstants.MAD_MEDIAN_RIGHT_OFFSET_KEY);
            //If the variant is near the read end (median threshold) and the positions are very similar (MAD threshold) then filter
            if ((tumorFwdPosMedian != null && tumorFwdPosMedian <= MTAC.PIR_MEDIAN_THRESHOLD && tumorFwdPosMAD != null && tumorFwdPosMAD <= MTAC.PIR_MAD_THRESHOLD) ||
                    (tumorRevPosMedian != null && tumorRevPosMedian <= MTAC.PIR_MEDIAN_THRESHOLD && tumorRevPosMAD != null && tumorRevPosMAD <= MTAC.PIR_MAD_THRESHOLD))
                filters.add(GATKVCFConstants.CLUSTERED_READ_POSITION_FILTER_NAME);
        }
    }

    private static void applyPanelOfNormalsFilter(final M2ArgumentCollection MTAC, final VariantContext vc, final Collection<String> filters) {
        final boolean siteInPoN = vc.hasAttribute(SomaticGenotypingEngine.IN_PON_VCF_ATTRIBUTE);
        if (siteInPoN) {
            filters.add(GATKVCFConstants.PON_FILTER_NAME);
        }
    }

    private static void applyGermlineVariantFilter(final M2ArgumentCollection MTAC, final VariantContext vc, final Collection<String> filters) {
        if (!vc.hasAttribute(GATKVCFConstants.NORMAL_LOD_KEY)) {
            return;
        }
        final boolean siteInCosmic = vc.hasAttribute(SomaticGenotypingEngine.IN_COSMIC_VCF_ATTRIBUTE);
        final boolean siteInDbsnp = vc.hasAttribute(SomaticGenotypingEngine.IN_DBSNP_VCF_ATTRIBUTE);
        if (siteInDbsnp && !siteInCosmic ) {
            // take the normal LOD of the first alt allele, which has the highest tumor LOD
            final double normalLod = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.NORMAL_LOD_KEY, () -> null, -1)[0];
            if (normalLod < MTAC.NORMAL_DBSNP_LOD_THRESHOLD) {
                filters.add(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME);
            }
        }
    }

    private static void applyStrandBiasFilter(final M2ArgumentCollection MTAC, final VariantContext vc, final Collection<String> filters) {
        if (MTAC.ENABLE_STRAND_ARTIFACT_FILTER) {
            if (vc.hasAttribute(GATKVCFConstants.TLOD_FWD_KEY) && vc.hasAttribute(GATKVCFConstants.TLOD_REV_KEY)
                    && vc.hasAttribute(GATKVCFConstants.TUMOR_SB_POWER_FWD_KEY) && vc.hasAttribute(GATKVCFConstants.TUMOR_SB_POWER_REV_KEY)) {
                final double forwardLod = vc.getAttributeAsDouble(GATKVCFConstants.TLOD_FWD_KEY, 0.0);
                final double reverseLod = vc.getAttributeAsDouble(GATKVCFConstants.TLOD_REV_KEY, 0.0);
                final double forwardPower = vc.getAttributeAsDouble(GATKVCFConstants.TUMOR_SB_POWER_FWD_KEY, 0.0);
                final double reversePower = vc.getAttributeAsDouble(GATKVCFConstants.TUMOR_SB_POWER_REV_KEY, 0.0);
                if ((forwardPower > MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && forwardLod < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD) ||
                        (reversePower > MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && reverseLod < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD)) {
                    filters.add(GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME);
                }
            }
        }
    }

    private static void applyEventDistanceFilters(final VariantContext vc, final Collection<String> filters) {
        final Integer eventCount = vc.getAttributeAsInt(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, -1);
        if (eventCount >= 3) {
            filters.add(GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME);
        }
    }

    //TODO: building a list via repeated side effects is ugly
    public static Set<String> calculateFilters(final M2ArgumentCollection MTAC, final VariantContext vc) {
        final Set<String> filters = new HashSet<>();
        applyEventDistanceFilters(vc, filters);
        applyTriallelicFilter(vc, filters);
        applyPanelOfNormalsFilter(MTAC, vc, filters);
        applyGermlineVariantFilter(MTAC, vc, filters);
        applyClusteredReadPositionFilter(MTAC, vc, filters);
        applyStrandBiasFilter(MTAC, vc, filters);
        applySTRFilter(vc, filters);

        return filters;
    }


}
