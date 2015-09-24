package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class CountReadsSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return CountReadsSpark.class.getSimpleName();
    }

    @Test
    public void test() throws Exception {
        final File unsortedBam = new File(getTestDataDir(), "count_reads.bam");
        final File outputTxt = createTempFile("count_reads", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(unsortedBam.getCanonicalPath());
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME); args.add(outputTxt.getCanonicalPath());
        this.runCommandLine(args.getArgsArray());

        final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile());
        Assert.assertEquals((int)Integer.valueOf(readIn), 8);
    }

}