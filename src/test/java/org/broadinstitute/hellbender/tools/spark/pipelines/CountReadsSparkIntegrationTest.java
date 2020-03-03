package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.charset.StandardCharsets;

public final class CountReadsSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return CountReadsSpark.class.getSimpleName();
    }


    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new Object[][]{
                {"count_reads.sam", null, 8L},
                {"count_reads.bam", null, 8L},
                {"count_reads.cram", "count_reads.fasta", 8L},
                {"count_reads_sorted.sam", null, 8L},
                {"count_reads_sorted.bam", null, 8L},
                {"count_reads_sorted.cram", "count_reads.fasta", 8L},
        };
    }

    @Test(dataProvider = "filenames", groups = "spark")
    public void testCountReads(final String fileIn, final String referenceName, final long expectedCount) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final File outputTxt = createTempFile("count_reads", ".txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(ORIG_BAM);
        args.addOutput(outputTxt);
        if (null != referenceName) {
            final File REF = new File(getTestDataDir(), referenceName);
            args.addReference(REF);
        }
        this.runCommandLine(args.getArgsArray());
        final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile(), StandardCharsets.UTF_8);
        Assert.assertEquals((long)Long.valueOf(readIn), expectedCount);
    }

    @Test(groups = "spark")
    public void test() throws Exception {
        final File unsortedBam = new File(getTestDataDir(), "count_reads.bam");
        final File outputTxt = createTempFile("count_reads", ".txt");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(unsortedBam);
        args.addOutput(outputTxt);
        this.runCommandLine(args.getArgsArray());

        final String readIn = FileUtils.readFileToString(outputTxt.getAbsoluteFile(), StandardCharsets.UTF_8);
        Assert.assertEquals((int)Integer.valueOf(readIn), 8);
    }

    @DataProvider(name="intervals")
    public Object[][] intervals(){
        return new Object[][]{
                new Object[]{"", 8l}, //no intervals specified see all reads including unaligned
                new Object[]{"-L chr7:1-20", 4l},
                new Object[]{"-L chr1", 0l},
                new Object[]{"-L chr1 -L chr7", 7l},
                new Object[]{"-XL chr7", 0l},
                new Object[]{"-XL chr7:2-404", 3l},
                new Object[]{"-L chr7:1-30 -L chr7:10-15 --" + IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME + " INTERSECTION", 3l},
                new Object[]{"-L chr7:1 --" + IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME + " 19", 4l },
                new Object[]{"-L " + getTestDataDir() + "/chr7_1_20_interval.list", 4l},
                new Object[]{"-L chr7:1-100 -XL chr7:2-100", 3l},
                new Object[]{"-L chr7:1-10 -L chr7:5-10 --"
                        +IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME  + " 10 --"
                        +IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME + " INTERSECTION -XL chr7:21-200", 4l }
        };
    }


    @Test(dataProvider = "intervals", groups = "spark")
    public void testCountReadsWithIntervals(final String interval_args, final long expectedCount) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), "count_reads_sorted.bam");
        final File outputFile = createTempFile("count_reads_spark","count");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(ORIG_BAM);
        args.addRaw(interval_args);
        args.addOutput(outputFile);

        this.runCommandLine(args.getArgsArray());

        try(XReadLines output = new XReadLines(outputFile)){
            Assert.assertEquals((long)Long.valueOf(output.next()), expectedCount);
        }
    }

    @Test(groups = "spark")
    public void testNoNPRWhenOutputIsUnspecified(){
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(new File(getTestDataDir(), "count_reads.bam"));
        this.runCommandLine(args.getArgsArray());
    }
}
