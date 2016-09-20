package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Set;

/**
 * VQSR-specific utility methods shared by VariantRecalibrator and ApplyVQSR.
 */
class VariantRecalibrationUtils {

    /**
     * Add the standard VCF header lines used with VQSR.
     * @param hInfo updated set of VCFHeaderLines
     */
    protected static void addVQSRStandardHeaderLines(final Set<VCFHeaderLine> hInfo) {
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VQS_LOD_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.CULPRIT_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSITIVE_LABEL_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NEGATIVE_LABEL_KEY));
    }

    /**
     * Add the standard allele-specific VCF header lines used with VQSR.
     * @param hInfo updated set of VCFHeaderLines
     */
    protected static void addAlleleSpecificVQSRHeaderLines(final Set<VCFHeaderLine> hInfo) {
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_FILTER_STATUS_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_CULPRIT_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VQS_LOD_KEY));
    }
}

