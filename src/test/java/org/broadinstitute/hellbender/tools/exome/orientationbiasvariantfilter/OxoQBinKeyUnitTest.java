package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SparkTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class OxoQBinKeyUnitTest extends BaseTest {

    @Test
    public void testCanSerializeOxoQBinKey() {

        final OxoQBinKey initialKey = OxoQScorer.createOxoQBinKey((byte) 'C', (byte) 'A', (byte) 'T', (byte) 'A', false);
        final OxoQBinKey k = SparkTestUtils.roundTripInKryo(initialKey, OxoQBinKey.class, SparkContextFactory.getTestSparkContext().getConf());
        Assert.assertEquals(k, initialKey);
    }

    @Test
    public void testOxoQBinKeyEquals() {

        final OxoQBinKey initialKey = OxoQScorer.createOxoQBinKey((byte) 'C', (byte) 'A', (byte) 'T', (byte) 'A', false);
        final OxoQBinKey k = OxoQScorer.createOxoQBinKey((byte) 'C', (byte) 'A', (byte) 'T', (byte) 'A', false);
        Assert.assertEquals(k, initialKey);
    }
}
