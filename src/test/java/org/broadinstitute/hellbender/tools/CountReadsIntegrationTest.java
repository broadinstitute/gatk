package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountReadsIntegrationTest extends CommandLineProgramTest {
    @Test(dataProvider = "filenames")
    public void testCountBases(String fileIn) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final String[] args = new String[]{
                "--input",  ORIG_BAM.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, 8l);
    }

    @DataProvider(name="filenames")
    public Object[][] filenames() {
        return new String[][]{
                {"count_reads.sam"},
                {"count_reads.bam"},
                {"count_reads_sorted.sam"},
                {"count_reads_sorted.bam"}
        };
    }

    @Test(dataProvider = "cram_filenames")
    public void testCountBasesCRAM(String fileIn, String refIn) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), fileIn);
        final File ref = new File(getTestDataDir(), refIn);
        final String[] args = new String[]{
                "--input",  ORIG_BAM.getAbsolutePath(),
                "--R", ref.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, 8l);
    }

    @DataProvider(name="cram_filenames")
    public Object[][] CRAMFilenames() {
        return new String[][]{
                {"count_reads_sorted.cram", "count_reads_sorted.fasta"},
                {"count_reads.cram", "count_reads.fasta"}
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
                new Object[]{"-L chr7:1-30 -L chr7:10-15 --interval_set_rule INTERSECTION", 3l},
                new Object[]{"-L chr7:1 --interval_padding 19", 4l },
                new Object[]{"-L " + getTestDataDir() + "/chr7_1_20.interval_list", 4l},
                new Object[]{"-L chr7:1-100 -XL chr7:2-100", 3l},
                new Object[]{"-L chr7:1-10 -L chr7:5-10 --interval_padding 10 --interval_set_rule INTERSECTION --XL chr7:21-200", 4l }
        };
    }


    @Test(dataProvider = "intervals")
    public void testCountBasesWithIntervals(String interval_args, long count) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), "count_reads_sorted.bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(ORIG_BAM.getAbsolutePath());
        args.add(interval_args);

        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, count);
    }

    // TODO: Fix and undisable the test cases that are disabled below.
    @DataProvider(name="intervalsWithCRAM")
    public Object[][] intervalsCRAM(){
        return new Object[][]{
                new Object[]{"-L 4:1", 0l},
                new Object[]{"", 11l}, //no intervals specified see all reads including unaligned
                new Object[]{"-L 4:1-20", 0l},
                new Object[]{"-L 1", 5l},
                new Object[]{"-L 1 -L 4", 6l},
                new Object[]{"-XL 4", 10l},
                new Object[]{"-XL 4:2-404", 11l},
                new Object[]{"-L 4:1-30 -L 4:10-15 --interval_set_rule INTERSECTION", 0l},
                new Object[]{"-L 4:1 --interval_padding 19", 0l },
                new Object[]{"-L " + getTestDataDir() + "/4_1_20.interval_list", 0l},
                new Object[]{"-L 4:1-100 -XL 4:2-100", 0l},
                new Object[]{"-L 4:1-10 -L 4:5-10 --interval_padding 10 --interval_set_rule INTERSECTION --XL 4:21-200", 0l }
        };
    }


    @Test(dataProvider = "intervalsWithCRAM")
    public void testCountBasesWithIntervalsCRAM(String interval_args, long count) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), "count_reads_sorted.cram");
        final File ref = new File(getTestDataDir(), "count_reads_sorted.fasta");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(ORIG_BAM.getAbsolutePath());
        args.add("--reference");
        args.add(ref.getAbsolutePath());
        args.add(interval_args);

        final Object res = this.runCommandLine(args.getArgsArray());
        Assert.assertEquals(res, count);
    }


}
