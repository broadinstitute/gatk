package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SVDeduplicatorUnitTest {

    final SVCollapser collapser = new SVCollapser(SVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END);
    final SVDeduplicator<SVCallRecord> deduplicator = new SVDeduplicator<>(collapser::collapse, SVTestUtils.dict);

    @Test
    public void testItemsAreIdentical() {
        //same bounds, different algs
        Assert.assertTrue(deduplicator.itemsAreIdentical(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff));

        //different bounds
        Assert.assertFalse(deduplicator.itemsAreIdentical(SVTestUtils.call1, SVTestUtils.call2));

        Assert.assertTrue(deduplicator.itemsAreIdentical(SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch));
    }

    @Test
    public void testDeduplicateItems() {
        final SVCallRecord merged1 = deduplicator.deduplicateSortedItems(Arrays.asList(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff)).get(0);
        Assert.assertEquals(merged1.getGenotypes().getSampleNames(), Sets.newHashSet(SVTestUtils.sample1.getSampleName(), SVTestUtils.sample2.getSampleName()));
        Assert.assertEquals(merged1.getAlgorithms().size(), 2);
        Assert.assertTrue(merged1.getAlgorithms().containsAll(SVTestUtils.depthAndStuff.getAlgorithms()));
    }

}