package org.broadinstitute.hellbender.tools.spark.bwa;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.pipelines.BwaAndMarkDuplicatesPipelineSpark;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class BwaAndMarkDuplicatesPipelineSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return BwaAndMarkDuplicatesPipelineSpark.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        //The expected results file was created by
        // 1) running BwaSpark on the input,
        // 2) running picard SortSam to sort by coordinate
        // 3) running picard MarkDuplicates on the result
        // 4) running picard SortSam to sort by queryname

        final File output = createTempFile("bwa", ".bam");
        if (!output.delete()) {
            Assert.fail();
        }

        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File ref = new File(b37_reference_20_21);
        args.addReference(ref);
        final File input = new File(largeFileTestDir, "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.queryname.noMD.bam");
        args.addInput(input);
        args.addOutput(output);
        this.runCommandLine(args.getArgsArray());

        final File expectedSam = new File(largeFileTestDir, "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.queryname.noMD.bwa.md.bam");
        SamAssertionUtils.assertSamsEqual(output, expectedSam);
    }

}
