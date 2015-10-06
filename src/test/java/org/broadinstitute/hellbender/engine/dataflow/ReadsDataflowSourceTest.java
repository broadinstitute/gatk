package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.Count;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public final class ReadsDataflowSourceTest extends BaseTest {
    // We also need GCS bucket and Spark tests (#666).

    private final File bam = new File(this.getToolTestDataDir(), "count_reads_sorted.bam");
    final String hiSeqBam = "src/test/resources/org/broadinstitute/hellbender/engine/dataflow/ReadsSource/HiSeq.1mb.1RG.2k_lines_invalid_sam_format.bam";

    @Test
    public void testGetHeaderFromLocalBAM(){
        final SAMFileHeader expectedHeader = SamReaderFactory.makeDefault().open(bam).getFileHeader();
        ReadsDataflowSource source = new ReadsDataflowSource(bam.getAbsolutePath(), null);
        SAMFileHeader header = source.getHeader();
        Assert.assertEquals(header, expectedHeader);
    }

    @Test
    public void testGetReadPCollectionLocal(){
        Pipeline p = GATKTestPipeline.create();
        ReadsDataflowSource source = new ReadsDataflowSource(bam.getAbsolutePath(), p);
        DataflowUtils.registerGATKCoders(p);
        PCollection<GATKRead> reads = source.getReadPCollection(ImmutableList.of(new SimpleInterval("chr7", 1, 404)), ValidationStringency.DEFAULT_STRINGENCY, false);
        PCollection<Long> count = reads.apply(Count.globally());
        DataflowAssert.thatSingleton(count).isEqualTo(7L);
        p.run();
    }

    // TODO: once ReadsSource has an option for including unmapped reads, we can enable this test.
    @Test(enabled = false)
    public void testLocalFile() {
        final String bam2 = "src/test/resources/org/broadinstitute/hellbender/tools/BQSR/HiSeq.1mb.1RG.2k_lines.alternate.bam";
        Pipeline pipeline = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(pipeline);
        ReadsDataflowSource readsSource = new ReadsDataflowSource(bam2, pipeline);
        SAMFileHeader header = readsSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        PCollection<GATKRead> reads = readsSource.getReadPCollection(intervals, ValidationStringency.SILENT, true);
        PCollection<Long> count = reads.apply(Count.globally());
        // for now we only get 1649, because it removes unmapped reads.
        DataflowAssert.thatSingleton(count).isEqualTo(1674L);
        pipeline.run();
    }


    //TODO renable this test issue #774
    @Test( enabled = false)
    public void testGetInvalidPCollectionLocal() {
        // ValidationStringency.SILENT should prevent any read error even though the input has what looks like invalid reads.
        Pipeline p = GATKTestPipeline.create();
        ReadsDataflowSource source = new ReadsDataflowSource(hiSeqBam, p);
        SAMFileHeader header = source.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        DataflowUtils.registerGATKCoders(p);
        PCollection<GATKRead> reads = source.getReadPCollection(IntervalUtils.getAllIntervalsForReference(sequenceDictionary), ValidationStringency.SILENT, false);
        PCollection<Long> count = reads.apply(Count.globally());
        // There are 1677 total reads in this file
        DataflowAssert.thatSingleton(count).isEqualTo(1677L);
        p.run();
    }
}
