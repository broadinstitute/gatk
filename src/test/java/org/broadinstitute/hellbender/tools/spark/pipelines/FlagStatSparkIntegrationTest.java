package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class FlagStatSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return FlagStatSpark.class.getSimpleName();
    }

    @Test
    public void flagStatSparkLocalNoInterval() throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--"+StandardArgumentDefinitions.INPUT_LONG_NAME); args.add( new File(getToolTestDataDir(),"flag_stat.bam").toString());
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        File outputFile = createTempFile("flagStatTest", ".txt");
        args.add(outputFile.getPath());

        this.runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists());
        //the expected output was created using stand-alone hellbender
        IntegrationTestSpec.assertMatchingFiles(Lists.newArrayList(outputFile), Lists.newArrayList(getToolTestDataDir() +"/"+ "expectedStats.txt"), false);
    }
    @Test
    public void flagStatSparkLocalWithBigInterval() throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();

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
        File outputFile = createTempFile("flagStatTest", ".txt");
        args.add(outputFile.getPath());

        this.runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists());
        //the expected output was created using stand-alone hellbender
        IntegrationTestSpec.assertMatchingFiles(Lists.newArrayList(outputFile), Lists.newArrayList(getToolTestDataDir() +"/"+ "expectedStats.chr1-chr8.txt"), false);
    }

    @Test
    public void flagStatSparkLocalWithSmallInterval() throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--"+StandardArgumentDefinitions.INPUT_LONG_NAME); args.add( new File(getToolTestDataDir(),"flag_stat.bam").toString());
        args.add("-L chr7:1-100 -XL chr7:2-100");
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        File outputFile = createTempFile("flagStatTest.chr1_1", ".txt");
        args.add(outputFile.getPath());

        this.runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists());
        //the expected output was created using stand-alone hellbender
        IntegrationTestSpec.assertMatchingFiles(Lists.newArrayList(outputFile), Lists.newArrayList(getToolTestDataDir() +"/"+ "expectedStats.chr1_1.txt"), false);
    }
}