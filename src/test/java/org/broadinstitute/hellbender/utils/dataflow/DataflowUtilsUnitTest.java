package org.broadinstitute.hellbender.utils.dataflow;

import com.cloudera.dataflow.spark.EvaluationResult;
import com.cloudera.dataflow.spark.SparkPipelineRunner;
import com.google.api.client.util.Lists;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.DoFnTester;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.ValidationStringency;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public final class DataflowUtilsUnitTest extends BaseTest {

    @Test
    public void testConvertToString(){
        Pipeline p = GATKTestPipeline.create();
        PCollection<Integer> pints = p.apply(Create.of(Arrays.asList(1, 2, 3)));
        PCollection<String> presults = pints.apply(DataflowUtils.convertToString());
        DataflowAssert.that(presults).containsInAnyOrder("1","2","3");
        p.run();
    }

    /**
     * Dataflow will wrap all exceptions into RuntimeExceptions
     */
    @Test(expectedExceptions = RuntimeException.class)
    public void testThrowException() {
        Pipeline p = GATKTestPipeline.create();
        PCollection<Integer> pints = p.apply(Create.of(Arrays.asList(1, 2, 3)));
        pints.apply(DataflowUtils.throwException(new IOException("fail")));
        p.run();
    }

    @Test
    public void testReadFromFileFn(){
        List<SimpleInterval> intervals = Arrays.asList(new SimpleInterval("chr7:1-202"), new SimpleInterval("chr8:2-202"));
        DoFn<File, GATKRead> readfn = new DataflowUtils.LoadReadsFromFileFn(intervals, ValidationStringency.SILENT);

        File inputFile = new File(getToolTestDataDir(), "example_reads.bam");
        List<GATKRead> expected = getReadsFromFile(intervals, inputFile, false);

        DoFnTester<File, GATKRead> tester = DoFnTester.of(readfn);
        List<GATKRead> output = tester.processBatch(inputFile);

        Assert.assertTrue(ReadUtils.readListsAreEqualIgnoreUUID(output, expected), "Actual reads do not match expected reads");
    }

    @Test
    public void testGetReadsFromHadoopBam() {
        List<SimpleInterval> intervals = Arrays.asList(new SimpleInterval("chr7:1-202"), new SimpleInterval("chr8:2-202"));
        File inputFile = new File(getToolTestDataDir(), "example_reads.bam");
        List<GATKRead> expected = getReadsFromFile(intervals, inputFile, true);

        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);
        PCollection<GATKRead> reads = DataflowUtils.getReadsFromHadoopBam(p, intervals, ValidationStringency.SILENT,
                new Path(inputFile.getAbsoluteFile().toURI()).toString());
        EvaluationResult result = SparkPipelineRunner.create().run(p);

        Assert.assertTrue(ReadUtils.readListsAreEqualIgnoreUUID(expected, Lists.newArrayList(result.get(reads))), "Actual reads do not match expected reads");
    }

    public List<GATKRead> getReadsFromFile(List<SimpleInterval> intervals, File inputFile, boolean useGoogleReads) {
        try(ReadsDataSource source = new ReadsDataSource(inputFile)) {
            source.setIntervalsForTraversal(intervals);

            return useGoogleReads ?
                    StreamSupport.stream(source.spliterator(), false)
                    .map(read -> new GoogleGenomicsReadToGATKReadAdapter(read.convertToGoogleGenomicsRead()))
                    .collect(Collectors.toList())
                    : StreamSupport.stream(source.spliterator(), false).collect(Collectors.toList());
        }
    }
}