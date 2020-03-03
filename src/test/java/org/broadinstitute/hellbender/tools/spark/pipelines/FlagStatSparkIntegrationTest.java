package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class FlagStatSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return FlagStatSpark.class.getSimpleName();
    }

    @Test(groups = "spark")
    public void flagStatSparkLocalNoInterval() throws IOException {
        File outputFile = createTempFile("flagStatTest", ".txt");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(getTestFile("flag_stat.bam"));
        args.addOutput(outputFile);

        this.runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists());
        //the expected output was created using stand-alone hellbender
        IntegrationTestSpec.assertMatchingFiles(Lists.newArrayList(outputFile), Lists.newArrayList(getToolTestDataDir() +"/"+ "expectedStats.txt"), false, null);
    }
    @Test(groups = "spark")
    public void flagStatSparkLocalWithBigInterval() throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addInput(getTestFile("flag_stat.bam"));
        args.addRaw("-L"); args.addRaw("chr1");
        args.addRaw("-L"); args.addRaw("chr2");
        args.addRaw("-L"); args.addRaw("chr3");
        args.addRaw("-L"); args.addRaw("chr4");
        args.addRaw("-L"); args.addRaw("chr5");
        args.addRaw("-L"); args.addRaw("chr6");
        args.addRaw("-L"); args.addRaw("chr7");
        args.addRaw("-L"); args.addRaw("chr8");
        File outputFile = createTempFile("flagStatTest", ".txt");
        args.addOutput(outputFile);

        this.runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists());
        //the expected output was created using stand-alone hellbender
        IntegrationTestSpec.assertMatchingFiles(Lists.newArrayList(outputFile), Lists.newArrayList(getToolTestDataDir() +"/"+ "expectedStats.chr1-chr8.txt"), false, null);
    }

    @Test(groups = "spark")
    public void flagStatSparkLocalWithSmallInterval() throws IOException {
        ArgumentsBuilder args = new ArgumentsBuilder();

        args.addInput(getTestFile("flag_stat.bam"));
        args.addRaw("-L chr7:1-100 -XL chr7:2-100");
        File outputFile = createTempFile("flagStatTest.chr1_1", ".txt");
        args.addOutput(outputFile);

        this.runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists());
        //the expected output was created using stand-alone hellbender
        IntegrationTestSpec.assertMatchingFiles(Lists.newArrayList(outputFile), Lists.newArrayList(getToolTestDataDir() +"/"+ "expectedStats.chr1_1.txt"), false, null);
    }

    @Test(groups = "spark")
    public void testNoNPRWhenOutputIsUnspecified(){
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(getTestFile("flag_stat.bam"));
        this.runCommandLine(args.getArgsArray());
    }
}
