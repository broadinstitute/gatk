package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.Count;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.BaseRecalOutput;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public final class ReadsSourceTest extends BaseTest {

    private final File bam = new File(this.getToolTestDataDir(), "count_reads_sorted.bam");
    final String hiSeqBam = "src/test/resources/org/broadinstitute/hellbender/engine/dataflow/ReadsSource/HiSeq.1mb.1RG.2k_lines.bam";

    @Test
    public void testGetHeaderFromLocalBAM(){
        final SAMFileHeader expectedHeader = SamReaderFactory.makeDefault().open(bam).getFileHeader();
        ReadsSource source = new ReadsSource(bam.getAbsolutePath(), null);
        SAMFileHeader header = source.getHeader();
        Assert.assertEquals(header, expectedHeader);
    }

    @Test
    public void testGetReadPCollectionLocal(){
        Pipeline p = GATKTestPipeline.create();
        ReadsSource source = new ReadsSource(bam.getAbsolutePath(), p);
        DataflowWorkarounds.registerGenomicsCoders(p);
        PCollection<Read> reads = source.getReadPCollection(ImmutableList.of(new SimpleInterval("chr7", 1, 404)), ValidationStringency.DEFAULT_STRINGENCY);
        PCollection<Long> count = reads.apply(Count.globally());
        DataflowAssert.thatSingleton(count).isEqualTo(7L);
        p.run();
    }
    /*
    // TODO: once ReadsSource has an option for including unmapped reads, we can add this test.
    @Test
    public void testLocalFile() {
        final String bam2 = "src/test/resources/org/broadinstitute/hellbender/tools/BQSR/HiSeq.1mb.1RG.2k_lines.alternate.bam";
        Pipeline pipeline = TestPipeline.create();
        DataflowWorkarounds.registerGenomicsCoders(pipeline);
        ReadsSource readsSource = new ReadsSource(bam2, pipeline);
        SAMFileHeader header = readsSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        PCollection<Read> reads = readsSource.getReadPCollection(intervals, ValidationStringency.SILENT);
        PCollection<Long> count = reads.apply(Count.globally());
        // for now we only get 1649, because it removes unmapped reads.
        DataflowAssert.thatSingleton(count).isEqualTo(1674L);
        pipeline.run();
    }
    */

    @Test
    public void testGetInvalidPCollectionLocal() {
        // ValidationStringency.SILENT should prevent any read error even though the input has what looks like invalid reads.
        Pipeline p = GATKTestPipeline.create();
        ReadsSource source = new ReadsSource(hiSeqBam, p);
        SAMFileHeader header = source.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        DataflowWorkarounds.registerGenomicsCoders(p);
        PCollection<Read> reads = source.getReadPCollection(IntervalUtils.getAllIntervalsForReference(sequenceDictionary), ValidationStringency.SILENT);
        PCollection<Long> count = reads.apply(Count.globally());
        DataflowAssert.thatSingleton(count).isEqualTo(1674L);
        p.run();
    }

    @Test(expectedExceptions = SAMException.class)
    public void testGetInvalidPCollectionLocalStrict() throws Throwable {
        // ValidationStringency.STRICT should trigger an error on an invalid file
        try {
            Pipeline p = GATKTestPipeline.create();
            ReadsSource source = new ReadsSource(hiSeqBam, p);
            SAMFileHeader header = source.getHeader();
            final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
            DataflowWorkarounds.registerGenomicsCoders(p);
            PCollection<Read> reads = source.getReadPCollection(IntervalUtils.getAllIntervalsForReference(sequenceDictionary), ValidationStringency.STRICT);
            PCollection<Long> count = reads.apply(Count.globally());
            p.run();
        } catch (RuntimeException x) {
            throw x.getCause();
        }
    }
}