package org.broadinstitute.hellbender.tools.dataflow;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


public class CountBasesFilteredIntegrationTest extends CommandLineProgramTest{



    @Test(groups = {"dataflow", "bucket"})
    public void testCountBasesFilteredIntegrationTestBucket(){

        List<String> args = new ArrayList<>();
        args.add("--bam"); args.add("gs://louisb_genomics/test/flag_stat.bam");
        args.add("--dataflowIntervals"); args.add("chr7:1-100000");
        args.add("--dataflowIntervals"); args.add("chr8:1-100000");
        args.add("--outputFile");
        File outputFile = createTempFile("countbasestest", ".txt");
        args.add(outputFile.getPath());

        runCommandLine(args);

        Assert.assertTrue(outputFile.exists());

    }

    @Test(groups = {"dataflow",})
    public void testCountBasesFilteredIntegrationTestLocal(){

        List<String> args = new ArrayList<>();
        args.add("--bam"); args.add(new File(getTestDataDir(),"flag_stat.bam").getPath());
        args.add("--dataflowIntervals"); args.add("chr7:1-100000");
        args.add("--dataflowIntervals"); args.add("chr8:1-100000");
        args.add("--outputFile");
        File outputFile = createTempFile("countbasestest", ".txt");
        args.add(outputFile.getPath());

        runCommandLine(args);

        Assert.assertTrue(outputFile.exists());

    }

    @Test(groups = {"dataflow","cloud","bucket"})
    public void testCountBasesFilteredInTheCloud(){
        List<String> args = new ArrayList<>();
        args.add("--bam"); args.add("gs://louisb_genomics/test/flag_stat.bam");
        args.add("--dataflowIntervals"); args.add("chr7:1-100000");
        args.add("--dataflowIntervals"); args.add("chr8:1-100000");
        args.add("--runner"); args.add("BLOCKING");
        args.add("--project"); args.add("broad-dsde-dev");
        args.add("--staging"); args.add("gs://louisb_genomics/staging");
        args.add("--outputFile");
        args.add("gs://louisb_genomics/test/output/testoutput.tmp.txt");

        runCommandLine(args);

    }

}