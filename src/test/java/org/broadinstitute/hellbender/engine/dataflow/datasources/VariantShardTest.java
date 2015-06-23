package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.UUID;

public final class VariantShardTest {

    @DataProvider(name = "variantShards")
    public Object[][] variantShards(){
        List<GATKRead> reads = Arrays.asList(makeRead(1, 300, 1), makeRead(100000, 10, 2), makeRead(299999, 2, 3));
        List<Iterable<VariantShard>> variantShards = Arrays.asList();
        variantShards.add(Arrays.asList(new VariantShard(0, "1")));
        variantShards.add(Arrays.asList(new VariantShard(1, "1")));
        variantShards.add(Arrays.asList(new VariantShard(2, "1"), new VariantShard(3, "1")));
        return new Object[][]{
                {reads, variantShards},
        };
    }

    public GATKRead makeRead(int start, int length, int i) {
        return ArtificialReadUtils.createRandomRead(start, length, i);
    }

    @Test(dataProvider = "variantShards")
    public void getVariantShardsFromIntervalTest(List<GATKRead> reads, List<Iterable<VariantShard>> shards) {
        for (int i = 0; i < reads.size(); ++i) {
            GATKRead r = reads.get(i);
            Iterable<VariantShard> expectedShards = shards.get(i);
            List<VariantShard> foundShards = VariantShard.getVariantShardsFromInterval(r);
            Assert.assertEquals(foundShards, expectedShards);
        }
    }
}