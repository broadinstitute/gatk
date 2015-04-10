package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Set;

//XXX this is a stub, placeholder to help other classes compile. May go away when we've put all the pieces together.
public final class VQSRUtils {
    public static void addVQSRStandardHeaderLines(final Set<VCFHeaderLine> hInfo) {
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VQS_LOD_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.CULPRIT_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSITIVE_LABEL_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NEGATIVE_LABEL_KEY));
    }
}
