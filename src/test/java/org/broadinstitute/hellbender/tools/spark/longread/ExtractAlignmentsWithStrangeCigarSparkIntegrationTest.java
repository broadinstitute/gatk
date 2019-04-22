package org.broadinstitute.hellbender.tools.spark.longread;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class ExtractAlignmentsWithStrangeCigarSparkIntegrationTest extends CommandLineProgramTest {

    @Test
    public void dummy() {
        final File inputBam = new File("/Users/shuang/Desktop/test.bam");//largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");
        final File outputBam = new File(gatkDirectory + "rr.bam");
        final List<String> argsList = Arrays.asList(new ArgumentsBuilder().addInput(inputBam).addOutput(outputBam).getArgsArray());
        runCommandLine(argsList);
    }
}
