package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;

public class SVConcordanceLinkageTest extends BaseTest {
    @Test
    public void testClusterTogetherVaryTypes() {
        final SVConcordanceLinkage linkage = new SVConcordanceLinkage(SVTestUtils.hg38Dict);
        for (final GATKSVVCFConstants.StructuralVariantAnnotationType type1 : GATKSVVCFConstants.StructuralVariantAnnotationType.values()) {
            // Pass in null strands to let them be determined automatically
            final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, SVTestUtils.getValidTestStrandA(type1),
                    "chr1", 2001, SVTestUtils.getValidTestStrandB(type1), type1, null, Collections.emptyList(),
                    SVTestUtils.getLength(1000, 2001, type1), Collections.emptyList(), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                    Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
            for (final GATKSVVCFConstants.StructuralVariantAnnotationType type2 : GATKSVVCFConstants.StructuralVariantAnnotationType.values()) {
                final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, SVTestUtils.getValidTestStrandA(type2),
                        "chr1", 2001, SVTestUtils.getValidTestStrandB(type2), type2, null, Collections.emptyList(),
                        SVTestUtils.getLength(1000, 2001, type2), Collections.emptyList(), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                // Should only cluster together if same type, except CNVs
                if ((type1 == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV && call2.isSimpleCNV()) ||
                        (type2 == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV && call1.isSimpleCNV())) {
                    Assert.assertTrue(linkage.areClusterable(call1, call2).getResult());
                } else {
                    Assert.assertEquals(linkage.areClusterable(call1, call2).getResult(), type1 == type2);
                }
            }
        }
    }
}