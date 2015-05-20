package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Count;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.common.collect.ImmutableList;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class ReadsSourceTest extends BaseTest {
    public static final String EXPECTED_HEADER = "@HD\tVN:1.0\tSO:coordinate\n" +
            "@SQ\tSN:chr1\tLN:101\n" +
            "@SQ\tSN:chr2\tLN:101\n" +
            "@SQ\tSN:chr3\tLN:101\n" +
            "@SQ\tSN:chr4\tLN:101\n" +
            "@SQ\tSN:chr5\tLN:101\n" +
            "@SQ\tSN:chr6\tLN:101\n" +
            "@SQ\tSN:chr7\tLN:404\n" +
            "@SQ\tSN:chr8\tLN:202\n" +
            "@RG\tID:0\tSM:Hi,Mom!\n" +
            "@PG\tID:1\tPN:Hey!\tVN:2.0\n";
    private final File bam = new File(this.getToolTestDataDir(), "count_reads_sorted.bam");

    @Test
    public void testGetHeaderFromLocalBAM(){
        ReadsSource source = new ReadsSource(bam.getAbsolutePath(), null);
        String header = source.getHeaderString();
        Assert.assertEquals(header, EXPECTED_HEADER);
    }

    @Test
    public void testGetReadPCollectionLocal(){
        Pipeline p = TestPipeline.create();
        ReadsSource source = new ReadsSource(bam.getAbsolutePath(), p);
        DataflowWorkarounds.registerGenomicsCoders(p);
        PCollection<Read> reads = source.getReadPCollection(ImmutableList.of(new SimpleInterval("chr7", 1, 404)));
        PCollection<Long> count = reads.apply(Count.globally());
        DataflowAssert.thatSingleton(count).isEqualTo(7L);
        p.run();
    }
}