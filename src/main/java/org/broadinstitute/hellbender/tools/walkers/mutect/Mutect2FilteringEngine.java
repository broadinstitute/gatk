package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collection;

/**
 * Created by davidben on 9/15/16.
 */
public class Mutect2FilteringEngine {
    final M2ArgumentCollection MTAC;

    public Mutect2FilteringEngine(final M2ArgumentCollection MTAC) {
        this.MTAC = MTAC;
    }

    public void applySTRFilter(final VariantContext vc, final Collection<String> filters) {
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

    public void applyClusteredReadPositionFilter(final VariantContext vc, final Collection<String> filters) {
        if (MTAC.ENABLE_CLUSTERED_READ_POSITION_FILTER) {
            final Double tumorFwdPosMedian = (Double) vc.getAttribute(Mutect2Engine.MEDIAN_LEFT_OFFSET_KEY);
            final Double tumorRevPosMedian = (Double) vc.getAttribute(Mutect2Engine.MEDIAN_RIGHT_OFFSET_KEY);
            final Double tumorFwdPosMAD = (Double) vc.getAttribute(Mutect2Engine.MAD_MEDIAN_LEFT_OFFSET_KEY);
            final Double tumorRevPosMAD = (Double) vc.getAttribute(Mutect2Engine.MAD_MEDIAN_RIGHT_OFFSET_KEY);
            //If the variant is near the read end (median threshold) and the positions are very similar (MAD threshold) then filter
            if ((tumorFwdPosMedian != null && tumorFwdPosMedian <= MTAC.PIR_MEDIAN_THRESHOLD && tumorFwdPosMAD != null && tumorFwdPosMAD <= MTAC.PIR_MAD_THRESHOLD) ||
                    (tumorRevPosMedian != null && tumorRevPosMedian <= MTAC.PIR_MEDIAN_THRESHOLD && tumorRevPosMAD != null && tumorRevPosMAD <= MTAC.PIR_MAD_THRESHOLD))
                filters.add(Mutect2Engine.CLUSTERED_READ_POSITION_FILTER_NAME);
        }
    }

    public void applyPanelOfNormalsFilter(final VariantContext vc, final Collection<String> filters, final FeatureContext featureContext) {
        //TODO: verify that vc.getStart() is the right second argument to use here
        final Collection<VariantContext> panelOfNormalsVC = featureContext.getValues(MTAC.normalPanelFeatureInput, vc.getStart());

        // Note that we're filtering even if the alleles are different.  This is due to the different roles
        // of the matched normal, which says which alternate *alleles* are germline events,
        // and the panel of normals, which says which *sites* are generally noisy and not to be trusted.
        if (!panelOfNormalsVC.isEmpty()) {
            filters.add(GATKVCFConstants.PON_FILTER_NAME);
        }
    }
}
