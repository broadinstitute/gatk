package org.broadinstitute.hellbender.tools.exome.conversion.allelicbalancecaller;

import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SparkTestUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class AllelicSplitCallerModelStateUnitTest extends BaseTest {

    @Test
    public void testBasicInit() {
        final ACNVModeledSegment acnvModeledSegment = new ACNVModeledSegment(new SimpleInterval("1", 1000, 1500),
                new PosteriorSummary(-4000, -4001, -4002),
                new PosteriorSummary(-4000, -4001, -4002));
        final List<ACNVModeledSegment> tempList = new ArrayList<>();
        tempList.add(acnvModeledSegment);

        final AllelicBalanceCallerModelState state = AllelicBalanceCallerModelState.createInitialCNLOHCallerModelState(0.2, tempList,
                HomoSapiensConstants.DEFAULT_PLOIDY, CNLOHCaller.NUM_RHOS);
        Assert.assertNotNull(state);
        Assert.assertNotNull(state.getEffectivePis());
        Assert.assertTrue(state.getEffectivePis().length > 0);
        Assert.assertTrue(state.getmVals().length > 0);
        Assert.assertTrue(state.getnVals().length > 0);
        Assert.assertEquals(MathUtils.sum(state.getEffectivePis()), 1.0, 1e-10);
    }

    @Test
    public void testSerializationRoundTrip() {
        final ACNVModeledSegment acnvModeledSegment = new ACNVModeledSegment(new SimpleInterval("1", 1000, 1500),
                new PosteriorSummary(-4000, -4001, -4002),
                new PosteriorSummary(-4000, -4001, -4002));
        final List<ACNVModeledSegment> tempList = new ArrayList<>();
        tempList.add(acnvModeledSegment);

        final AllelicBalanceCallerModelState state = AllelicBalanceCallerModelState.createInitialCNLOHCallerModelState(0.2, tempList,
                HomoSapiensConstants.DEFAULT_PLOIDY, CNLOHCaller.NUM_RHOS);
        SparkTestUtils.roundTripInKryo(state, AllelicBalanceCallerModelState.class,
                SparkContextFactory.getTestSparkContext().getConf());
    }
}
