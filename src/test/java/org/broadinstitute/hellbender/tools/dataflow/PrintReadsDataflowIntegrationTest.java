package org.broadinstitute.hellbender.tools.dataflow;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.testng.Assert.*;

public class PrintReadsDataflowIntegrationTest extends CommandLineProgramTest {
    @Test(groups = {"dataflow", "bucket"})
    public void testPrintReadsDataflowIntegrationTest() throws IOException {

        List<String> args = new ArrayList<>();
        args.add("--bam"); args.add("gs://louisb_genomics/test/flag_stat.bam");
        args.add("--dataflowIntervals"); args.add("chr7:1-100000");
        args.add("--dataflowIntervals"); args.add("chr8:1-100000");
        args.add("--outputFile");
        File outputFile = createTempFile("printreadstest", ".txt");
        args.add(outputFile.getPath());

        runCommandLine(args);


     //   IntegrationTestSpec.assertMatchingFiles(outputFile, );

    }

    @Test
    public void testPackage(){
        System.err.println(this.getClass().getPackage().toString());
        Assert.fail();
    }
}