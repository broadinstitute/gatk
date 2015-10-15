package org.broadinstitute.hellbender.tools.spark.validation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class CompareMDIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test() throws Exception {
        //final File firstBam = new File("/Users/davidada/dev/git_projects/broad/hellbender/deduped-CEUTrio.HiSeq.WGS.b37.ch20.4m-6m.NA12878.bam");
        //final File secondBam = new File("/Users/davidada/dev/git_projects/broad/hellbender/gatk-deduped-CEUTrio.HiSeq.WGS.b37.ch20.4m-6m.NA12878.bam/");
        //final File secondBam = new File("/Users/davidada/dev/git_projects/broad/hellbender/gatk-deduped-tmp-mixed-CEUTrio.HiSeq.WGS.b37.ch20.4m-6m.NA12878.bam/");
        //final File secondBam = new File("/Users/davidada/dev/git_projects/broad/hellbender/gatk-deduped-tmp-mixed2-CEUTrio.HiSeq.WGS.b37.ch20.4m-6m.NA12878.bam/");

        final File firstBam = new File("/Users/davidada/dev/git_projects/broad/hellbender/deduped-CEUTrio.HiSeq.WGS.b37.ch20.4m-12m.NA12878.bam");
        final File secondBam = new File("/Users/davidada/dev/git_projects/broad/hellbender/gatk-deduped-tmp-mixed2-CEUTrio.HiSeq.WGS.b37.ch20.4m-12m.NA12878.bam");

        //

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(firstBam.getCanonicalPath());
        args.add("--" + "I2");
        args.add(secondBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());
        Assert.assertTrue(false);
    }

}