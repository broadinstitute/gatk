package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
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

    @DataProvider(name="intervals")
    public Object[][] intervals(){
        return new Object[][]{
                new Object[]{"",8L}, // no intervals specified, see all reads that are aligned and not aligned
                new Object[]{"-L chr7:1", 3l},
                new Object[]{"-L chr7:1-20", 4l},
                new Object[]{"-L chr1", 0l},
                new Object[]{"-L chr1 -L chr7", 7l},
                new Object[]{"-XL chr7", 0l},
                new Object[]{"-XL chr7:2-404", 3l},
                new Object[]{"-L chr7:1-30 -L chr7:10-15 --interval_set_rule INTERSECTION", 3l},
                new Object[]{"-L chr7:1 --interval_padding 19", 4l },
                new Object[]{"-L " + getTestDataDir() + "/chr7_1_20.interval_list", 4l},
                new Object[]{"-L chr7:1-100 -XL chr7:2-100", 3l},
                new Object[]{"-L chr7:1-10 -L chr7:5-10 --interval_padding 10 --interval_set_rule INTERSECTION --XL chr7:21-200", 4l }
        };
    }


    @Test(dataProvider = "intervals")
    public void testCountReadsWithIntervals(final String interval_args, final long expectedCount) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), "count_reads_sorted.bam");
        final File outputFile = createTempFile("count_reads_spark","count");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(ORIG_BAM.getAbsolutePath());
        args.add(interval_args);
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME); args.add(outputFile);

        this.runCommandLine(args.getArgsArray());

        try(XReadLines output = new XReadLines(outputFile)){
            Assert.assertEquals((long)Long.valueOf(output.next()), expectedCount);
        }

    }
}