package org.broadinstitute.hellbender.tools.spark.bwa;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class BwaSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return BwaSpark.class.getSimpleName();
    }

    @Test
    public void testPairedEnd() throws Exception {
        final File expectedSam = getTestFile("bwa.sam");

        final File ref = getTestFile("ref.fa");
        final File input = getTestFile("R.bam"); // this is a queryname sorted, paired, and unaligned input

        final File output = createTempFile("bwa", ".bam");
        Assert.assertTrue(output.delete());

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, ref);
        args.addFileArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, input);
        args.addBooleanArgument(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME, true); // disable since input does not have a sequence dictionary
        args.addBooleanArgument(GATKSparkTool.SHARDED_OUTPUT_LONG_NAME, true);
        args.addArgument(GATKSparkTool.NUM_REDUCERS_LONG_NAME,"1");
        args.addOutput(output);
        args.addFileArgument(BwaArgumentCollection.BWA_MEM_INDEX_IMAGE_FULL_NAME, getTestFile("ref.fa.img"));
        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(new File(output, "part-r-00000.bam"), expectedSam);
    }

    @Test
    public void testSingleEnd() throws Exception {
        final File expectedSam = getTestFile("seBwa.bam");

        final File ref = getTestFile("ref.fa");
        final File input = getTestFile("seR.bam");
        final File output = createTempFile("bwa", ".bam");
        Assert.assertTrue(output.delete());

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, ref);
        args.addFileArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, input);
        args.addBooleanArgument(StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME , true); // disable since input does not have a sequence dictionary
        args.addBooleanArgument(GATKSparkTool.SHARDED_OUTPUT_LONG_NAME, true);
        args.addArgument(GATKSparkTool.NUM_REDUCERS_LONG_NAME,"1");
        args.addOutput(output);
        args.add("--" + BwaArgumentCollection.SINGLE_END_ALIGNMENT_FULL_NAME);
        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(new File(output, "part-r-00000.bam"), expectedSam);
    }
}
