package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Count;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class ReadsSourceTest extends BaseTest {

    private final File bam = new File(this.getToolTestDataDir(), "count_reads_sorted.bam");

    @Test
    public void testGetHeaderFromLocalBAM(){
        final SAMFileHeader expectedHeader = SamReaderFactory.makeDefault().open(bam).getFileHeader();
        ReadsSource source = new ReadsSource(bam.getAbsolutePath(), null);
        SAMFileHeader header = source.getHeader();
        Assert.assertEquals(header, expectedHeader);
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