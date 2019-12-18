package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SVClusterEngineTest {

    private final SVClusterEngine engine = new SVClusterEngine(SVTestUtils.dict);

    @BeforeTest
    public void initializeDefragmenters() {
        engine.add(SVTestUtils.call1);
    }

    @Test
    public void testFlattenCluster() {
        //depth only and depthAndStuff have same bounds, less than call2
        final List<SVCallRecordWithEvidence> testCluster = Arrays.asList(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff, SVTestUtils.call2);
        final SVCallRecordWithEvidence flattened = engine.flattenCluster(testCluster);
        Assert.assertEquals(flattened.getStart(), SVTestUtils.depthAndStuff.getStart());
        Assert.assertEquals(flattened.getEnd(), SVTestUtils.depthAndStuff.getEnd());
        //should have all the algs
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.depthAndStuff.getAlgorithms()));
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.depthOnly.getAlgorithms()));
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.call2.getAlgorithms()));
        //should have all the genotypes
        Assert.assertTrue(flattened.getGenotypes().containsAll(SVTestUtils.depthAndStuff.getGenotypes()));
        Assert.assertTrue(flattened.getGenotypes().containsAll(SVTestUtils.depthOnly.getGenotypes()));
        Assert.assertTrue(flattened.getGenotypes().containsAll(SVTestUtils.call2.getGenotypes()));
        //TODO: add test for insertion cluster
    }

    @Test
    public void testClusterTogether() {
        Assert.assertTrue(engine.clusterTogether(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff));
        Assert.assertFalse(engine.clusterTogether(SVTestUtils.depthOnly, SVTestUtils.inversion));
        Assert.assertFalse(engine.clusterTogether(SVTestUtils.call1, SVTestUtils.call2));
    }

    @Test
    public void testGetClusteringInterval() {
        Assert.assertTrue(engine.getClusteringInterval(SVTestUtils.leftEdgeCall, null).getStart() > 0);
        Assert.assertTrue(engine.getClusteringInterval(SVTestUtils.rightEdgeCall, null).getEnd() < SVTestUtils.chr1Length);

        final SimpleInterval littleCluster = engine.getClusteringInterval(SVTestUtils.call1, null);
        final SimpleInterval totalInterval = engine.getClusteringInterval(SVTestUtils.call2, littleCluster);
        //TODO: add more quantitative checks
        //min start for combined interval should be greater than the leftmost bound, which is the start of call1
        Assert.assertTrue(totalInterval.getStart() < SVTestUtils.call1.getStart());
        //max start for combined interval should be greater than the leftmost bound, and less than the nearest event end
        Assert.assertTrue(totalInterval.getEnd() > SVTestUtils.call1.getStart());
        Assert.assertTrue(totalInterval.getEnd() > SVTestUtils.call1.getEnd());
    }

    @Test
    public void testItemsAreIdentical() {
        //same bounds, different algs
        Assert.assertTrue(engine.itemsAreIdentical(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff));

        //different bounds
        Assert.assertFalse(engine.itemsAreIdentical(SVTestUtils.call1, SVTestUtils.call2));
    }

    @Test
    public void testDeduplicateIdenticalItems() {
        final SVCallRecordWithEvidence merged1 = engine.deduplicateIdenticalItems(Arrays.asList(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff));
        Assert.assertEquals(merged1.getGenotypes().size(), 2);
        Assert.assertTrue(merged1.getGenotypes().containsAll(Arrays.asList(SVTestUtils.sample1, SVTestUtils.sample2)));
        Assert.assertEquals(merged1.getAlgorithms().size(), 2);
        Assert.assertTrue(merged1.getAlgorithms().containsAll(SVTestUtils.depthAndStuff.getAlgorithms()));
    }

    @Test
    public void testIsDepthOnlyCall() {
        Assert.assertTrue(SVClusterEngine.isDepthOnlyCall(SVTestUtils.depthOnly));
        Assert.assertFalse(SVClusterEngine.isDepthOnlyCall(SVTestUtils.depthAndStuff));
    }
}