package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.appengine.repackaged.com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class PrintReadsDataflowIntegrationTest extends CommandLineProgramTest {

    @Test(groups = {"dataflow"})
    public void testPrintReadsDataflowIntegrationTest() throws IOException {

        List<String> args = new ArrayList<>();
        File inputFile = new File(getToolTestDataDir()+"/"+ "print_reads.bam");
        File placeHolder = createTempFile("printreadstest", ".txt");
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(inputFile.toString());
        args.add("--L"); args.add("chr7:1-202");
        args.add("--L"); args.add("chr8:1-202");
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME); args.add(placeHolder.getPath());

        runCommandLine(args);
        File outputFile = findDataflowOutput(placeHolder);

        Set<String> expectedReadStrings;

        try(ReadsDataSource readsExpected = new ReadsDataSource(inputFile)){
            readsExpected.setIntervalsForTraversal(Lists.newArrayList(new SimpleInterval("chr7:1-202"), new SimpleInterval("chr8:1-202")));
            expectedReadStrings = StreamSupport.stream(readsExpected.spliterator(), true)
                    .map(SAMRecord::getSAMString)
                    .map(String::trim)
                    .collect(Collectors.toSet());
        }

        Set<String> foundReadStrings;
        try(XReadLines outputStrings = new XReadLines(outputFile)){
            foundReadStrings = Sets.newHashSet(outputStrings.iterator());
        }

        Assert.assertEquals(foundReadStrings, expectedReadStrings);

    }

}