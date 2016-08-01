package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class CalculateCoverageSparkIntegrationTest extends CommandLineProgramTest {

    @Test(groups = "spark")
    public void test() throws Exception {
        final File bam = new File(publicTestDir + "large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");
        final File out = File.createTempFile("out", ".txt");
        out.delete();
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(bam.getAbsolutePath());
        args.add("--output");
        args.add(out.getAbsolutePath());
        this.runCommandLine(args.getArgsArray());

        System.out.println("Output: " + out);
    }
}
