package org.broadinstitute.hellbender.utils.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.appengine.repackaged.com.google.common.collect.Lists;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.DoFnTester;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class DataflowUtilsUnitTest extends BaseTest {

    @Test
    public void testConvertToString(){
        Pipeline p = TestPipeline.create();
        PCollection<Integer> pints = p.apply(Create.of(Lists.newArrayList(1, 2, 3)));
        PCollection<String> presults = pints.apply(DataflowUtils.convertToString());
        DataflowAssert.that(presults).containsInAnyOrder("1","2","3");
        p.run();
    }

    /**
     * Dataflow will wrap all exceptions into RuntimeExceptions
     */
    @Test(expectedExceptions = RuntimeException.class)
    public void testThrowException() {
        Pipeline p = TestPipeline.create();
        PCollection<Integer> pints = p.apply(Create.of(Lists.newArrayList(1, 2, 3)));
        pints.apply(DataflowUtils.throwException(new IOException("fail")));
        p.run();
    }

    @Test
    public void testReadFromFileFn(){
        List<SimpleInterval> intervals = Lists.newArrayList(new SimpleInterval("chr7:1-202"), new SimpleInterval("chr8:2-202"));
        DoFn<File, Read> readfn = new DataflowUtils.LoadReadsFromFileFn(intervals);

        File inputFile = new File(getToolTestDataDir(), "example_reads.bam");
        List<Read> expected = getReadsFromFile(intervals, inputFile);

        DoFnTester<File, Read> tester = DoFnTester.of(readfn);
        List<Read> output = tester.processBatch(inputFile);

        Assert.assertEquals(output, expected);
    }

    public List<Read> getReadsFromFile(List<SimpleInterval> intervals, File inputFile) {
        try(ReadsDataSource source = new ReadsDataSource(inputFile)) {
            source.setIntervalsForTraversal(intervals);

            return StreamSupport.stream(source.spliterator(), false)
                    .map(ReadConverter::makeRead)
                    .collect(Collectors.toList());
        }
    }
}