package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class ReferenceShardUnitTest extends GATKBaseTest {

    @DataProvider(name = "reads")
    public Object[][] reads() {
        return new Object[][]{
                {ReadsPreprocessingPipelineTestData.makeRead("1", 1, 300, 1, SAMRecord.class), new ReferenceShard(0, "1")},  //right in the middle of the shard
                {ReadsPreprocessingPipelineTestData.makeRead("1", ReferenceShard.REFERENCE_SHARD_SIZE, 10, 2, SAMRecord.class), new ReferenceShard(1, "1")}, //at  the start of a shard
                {ReadsPreprocessingPipelineTestData.makeRead("1", 3 * ReferenceShard.REFERENCE_SHARD_SIZE - 1, 2, 3, SAMRecord.class),  new ReferenceShard(2, "1")}  //overlapping the end of a shard
        };
    }

    @Test(dataProvider = "reads")
    public void getVariantShardsFromIntervalTest(GATKRead r, ReferenceShard expectedShard) {
        ReferenceShard foundShard = ReferenceShard.getShardNumberFromInterval(r);
        Assert.assertEquals(foundShard, expectedShard);
    }
}