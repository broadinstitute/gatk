package org.broadinstitute.hellbender.tools.exome.cnlohcaller;

import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.SerializationTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class CNLOHCallerModelStateUnitTest extends BaseTest {

    @Test
    public void testBasicInit() {
        final ACNVModeledSegment acnvModeledSegment = new ACNVModeledSegment(new SimpleInterval("1", 1000, 1500),
                new PosteriorSummary(-4000, -4001, -4002),
                new PosteriorSummary(-4000, -4001, -4002));
        final List<ACNVModeledSegment> tempList = new ArrayList<>();
        tempList.add(acnvModeledSegment);

        final CNLOHCallerModelState state = CNLOHCallerModelState.createInitialCNLOHCallerModelState(0.2, tempList);
        Assert.assertNotNull(state);
        Assert.assertNotNull(state.getCnToPiMap());
        Assert.assertTrue(state.getCnToPiMap().entrySet().size() > 0);

        final List<Integer> keysAsList = new ArrayList<>(state.getCnToPiMap().keySet());
        Collections.sort(keysAsList);
        Assert.assertEquals(state.getCnToPiMap().keySet(), keysAsList);
    }

    @Test
    public void testSerializationRoundTrip() {
        final ACNVModeledSegment acnvModeledSegment = new ACNVModeledSegment(new SimpleInterval("1", 1000, 1500),
                new PosteriorSummary(-4000, -4001, -4002),
                new PosteriorSummary(-4000, -4001, -4002));
        final List<ACNVModeledSegment> tempList = new ArrayList<>();
        tempList.add(acnvModeledSegment);

        final CNLOHCallerModelState state = CNLOHCallerModelState.createInitialCNLOHCallerModelState(0.2, tempList);
        SerializationTestUtils.roundTripInKryo(state, CNLOHCallerModelState.class,
                SparkContextFactory.getTestSparkContext().getConf());
    }
}
