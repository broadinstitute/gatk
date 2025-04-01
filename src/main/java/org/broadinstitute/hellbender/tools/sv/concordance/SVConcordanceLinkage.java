package org.broadinstitute.hellbender.tools.sv.concordance;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;

public class SVConcordanceLinkage extends CanonicalSVLinkage<SVCallRecord> {

    public SVConcordanceLinkage(final SAMSequenceDictionary dictionary) {
        super(dictionary, false);
    }

    @Override
    protected boolean cnvTypesMatch(final SVCallRecord a, final SVCallRecord b) {
        // Allow multi-allelic CNVs to cluster with both DELs and DUPs
        return (a.isSimpleCNV() && b.isSimpleCNV()) &&
                (a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV ||
                        b.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV);
    }
}
