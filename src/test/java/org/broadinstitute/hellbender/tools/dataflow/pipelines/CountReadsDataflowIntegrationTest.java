package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.common.base.Strings;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CountReadsDataflowIntegrationTest extends DataflowCommandLineProgramTest {

    @DataProvider(name="intervals")
    public Object[][] intervals(){
        return new Object[][]{
                new Object[]{"",7L}, // no intervals specified, see all reads that are aligned
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
    public void testCountBasesWithIntervals(String interval_args, long count) throws Exception {
        final File ORIG_BAM = new File(getTestDataDir(), "count_reads_sorted.bam");
        final File placeHolder = createTempFile("count_reads_dataflow","count");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(ORIG_BAM.getAbsolutePath());
        args.add(interval_args);
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME); args.add(placeHolder);
        addDataflowRunnerArgs(args);

        this.runCommandLine(args.getArgsArray());

        File outputFile = findDataflowOutput(placeHolder);
        try(XReadLines output = new XReadLines(outputFile)){
            Assert.assertEquals((long)Long.valueOf(output.next()), count);
        }

    }

}