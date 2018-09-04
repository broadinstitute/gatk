package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountReadsIntegrationTest extends CommandLineProgramTest {

    @Test(dataProvider = "filenames")
    public void testCountReads(final String fileIn, final String referenceName) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(ORIG_BAM.getAbsolutePath());
        if (null != referenceName) {
            final File REF = new File(getTestDataDir(), referenceName);
            args.add("-R");
            args.add(REF.getAbsolutePath());
        }
        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, 8l);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"count_reads.sam", null},
                {"count_reads.bam", null},
                {"count_reads.cram", "count_reads.fasta"},
                {"count_reads_sorted.sam", null},
                {"count_reads_sorted.bam", null},
                {"count_reads_sorted.cram", "count_reads.fasta"},
        };
    }

    @DataProvider(name="intervals")
    public Object[][] intervals(){
        return new Object[][]{
                new Object[]{"-L chr7:1", 3l},
                new Object[]{"", 8l}, //no intervals specified see all reads including unaligned
                new Object[]{"-L chr7:1-20", 4l},
                new Object[]{"-L chr1", 0l},
                new Object[]{"-L chr1 -L chr7", 7l},
                new Object[]{"-XL chr7", 0l},
                new Object[]{"-XL chr7:2-404", 3l},
                new Object[]{"-L chr7:1-30 -L chr7:10-15 --" + IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME + " INTERSECTION", 3l},
                new Object[]{"-L chr7:1 --" + IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME + " 19", 4l },
                new Object[]{"-L " + getTestDataDir() + "/chr7_1_20.interval_list", 4l},
                new Object[]{"-L chr7:1-100 -XL chr7:2-100", 3l},
                new Object[]{"-L chr7:1-10 -L chr7:5-10 --"
                        +IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME  + " 10 --"
                        +IntervalArgumentCollection.INTERVAL_SET_RULE_LONG_NAME + " INTERSECTION -XL chr7:21-200", 4l }
        };
    }

    @Test(dataProvider = "intervals")
    public void testCountBAMReadsWithIntervals(final String interval_args, final long count) throws Exception {
        countReads(interval_args, "count_reads_sorted.bam", null, count);
    }

    @Test(dataProvider = "intervals")
    public void testCountCRAMReadsWithIntervals(final String interval_args, final long count) throws Exception {
        countReads(interval_args, "count_reads_sorted.cram", "count_reads.fasta", count);
    }

    private void countReads(final String interval_args, final String fileName, final String referenceName, final long count) {
        final File ORIG_BAM = new File(getTestDataDir(), fileName);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(ORIG_BAM.getAbsolutePath());
        if (null != referenceName) {
            final File REF = new File(getTestDataDir(), referenceName);
            args.add("-R");
            args.add(REF.getAbsolutePath());
        }
        args.add(interval_args);

        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, count);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testNoArgs() throws Exception {
        //making sure that this blows up in a very specific way (UserException.CommandLineException) if we give bogus arguments
        final ArgumentsBuilder args = new ArgumentsBuilder();
        this.runCommandLine(args.getArgsArray());
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testBogusArgs() throws Exception {
        //making sure that this blows up in a very specific way (UserException.CommandLineException) if we give bogus arguments
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--fred");
        this.runCommandLine(args.getArgsArray());
    }
}
