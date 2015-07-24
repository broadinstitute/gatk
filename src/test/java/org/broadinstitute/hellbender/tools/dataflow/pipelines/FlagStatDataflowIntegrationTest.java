package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public final class FlagStatDataflowIntegrationTest extends CommandLineProgramTest {


    @Test( groups = {"bucket", "dataflow"})
    public void flagStatDataflowLocalReadFromBucket() throws IOException {
        List<String> args = new ArrayList<>();
        args.add("--apiKey"); args.add(getDataflowTestApiKey());
        args.add("--project"); args.add(getDataflowTestProject());
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME); args.add( getDataflowTestInputPath() + "org/broadinstitute/hellbender/tools/dataflow/FlagStatDataflow/flag_stat.bam");
        args.add("-L"); args.add("chr1");
        args.add("-L"); args.add("chr2");
        args.add("-L"); args.add("chr3");
        args.add("-L"); args.add("chr4");
        args.add("-L"); args.add("chr5");
        args.add("-L"); args.add("chr6");
        args.add("-L"); args.add("chr7");
        args.add("-L"); args.add("chr8");
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        File placeHolder = createTempFile("flagStatTest", ".txt");
        args.add(placeHolder.getPath());

        addDataflowRunnerArgs(args);


        runCommandLine(args);

        File outputFile = findDataflowOutput(placeHolder);
        Assert.assertTrue(outputFile.exists());
        IntegrationTestSpec.assertMatchingFiles(Lists.newArrayList(outputFile), Lists.newArrayList(getToolTestDataDir()+"/" + "expectedStats.txt"));


    }

    @Test(groups = "dataflow")
    public void flagStatDataflowLocal() throws IOException {
        List<String> args = new ArrayList<>();
        args.add("--"+StandardArgumentDefinitions.INPUT_LONG_NAME); args.add( new File(getToolTestDataDir(),"flag_stat.bam").toString());
        args.add("-L"); args.add("chr1");
        args.add("-L"); args.add("chr2");
        args.add("-L"); args.add("chr3");
        args.add("-L"); args.add("chr4");
        args.add("-L"); args.add("chr5");
        args.add("-L"); args.add("chr6");
        args.add("-L"); args.add("chr7");
        args.add("-L"); args.add("chr8");
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        File placeHolder = createTempFile("flagStatTest", ".txt");
        args.add(placeHolder.getPath());

        addDataflowRunnerArgs(args);


        runCommandLine(args);

        File outputFile = findDataflowOutput(placeHolder);
        Assert.assertTrue(outputFile.exists());
        IntegrationTestSpec.assertMatchingFiles(Lists.newArrayList(outputFile), Lists.newArrayList(getToolTestDataDir() +"/"+ "expectedStats.txt"));


    }



}