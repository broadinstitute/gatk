package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

import static org.testng.Assert.*;

public class SVClusterEngineTest {

    private final SVCallRecordWithEvidence depthOnly = new SVCallRecordWithEvidence("chr1", 10000, true, "chr1", 20000, true,
            StructuralVariantType.CNV, 10001, Arrays.asList("depth"), Collections.emptyList(),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    private final SVCallRecordWithEvidence depthAndStuff = new SVCallRecordWithEvidence("chr1", 10000, true, "chr1", 20000, true,
            StructuralVariantType.CNV, 10001, Arrays.asList("depth", "PE"), Collections.emptyList(),
            Collections.emptyList(), Collections.emptyList(), Collections.emptyList());

    @Test
    public void testFlattenCluster() {
    }

    @Test
    public void testClusterTogether() {
    }

    @Test
    public void testGetClusteringInterval() {
    }

    @Test
    public void testItemsAreIdentical() {
    }

    @Test
    public void testDeduplicateIdenticalItems() {
    }

    @Test
    public void testIsDepthOnlyCall() {
    }
}