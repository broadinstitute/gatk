package org.broadinstitute.hellbender.tools.spark.bwa;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.pipelines.BwaAndMarkDuplicatesPipelineSpark;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class BwaAndMarkDuplicatesPipelineSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return BwaAndMarkDuplicatesPipelineSpark.class.getSimpleName();
    }

    @Test
    public void test() throws Exception {
        //This file was created by 1) running bwaspark on the input and 2) running picard MarkDuplicates on the result
        final File expectedSam = new File(largeFileTestDir, "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.md.bam");

        final File ref = new File(b37_reference_20_21);
        final File input = new File(largeFileTestDir, "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.unaligned.bam");
        final File output = createTempFile("bwa", ".bam");
        if (!output.delete()) {
            Assert.fail();
        }

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(ref);
        args.addInput(input);
        args.addOutput(output);
        args.addBooleanArgument(MarkDuplicatesSparkArgumentCollection.DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME, true);
        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(output, expectedSam);
    }

}
