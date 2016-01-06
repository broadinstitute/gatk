package org.broadinstitute.hellbender.tools.spark.validation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class CompareDuplicatesSparkIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name = "CompareIdenticalDuplicatesProvider")
    public Object[][] makeCompareIdenticalDuplicatesProvider() {
        final String resourceDir = getTestDataDir() + "/validation/";
        final File legacyBam = new File(resourceDir, "tmp.legacy.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");
        final File gatk4Bam = new File(resourceDir, "tmp.gatk4.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");
        return new Object[][]{
                {legacyBam, legacyBam},
                {legacyBam, gatk4Bam},
        };
    }

    @DataProvider(name = "CompareDifferentDuplicatesProvider")
    public Object[][] makeCompareDifferentDuplicatesProvider() {
        final String resourceDir = getTestDataDir() + "/validation/";
        final File gatk4Bam = new File(resourceDir, "tmp.gatk4.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");
        final File singleReadBam = new File(resourceDir, "single.read.bam");
        final File ceuTrioBam = new File(resourceDir, "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");
        final File dupfreeBam = new File(resourceDir, "clean.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");

        return new Object[][]{
                {gatk4Bam, singleReadBam},
                {gatk4Bam, ceuTrioBam},
                {gatk4Bam, dupfreeBam},
        };
    }

    @Test(dataProvider = "CompareIdenticalDuplicatesProvider")
    public void identicalBamTest(File firstBam, File secondBam) throws Exception {
        // These files are the same and should produce no diffs.

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(firstBam.getCanonicalPath());
        args.add("--" + "I2");
        args.add(secondBam.getCanonicalPath());
        args.add("--" + "throwOnDiff");
        args.add("true");

        this.runCommandLine(args.getArgsArray());
    }

    @Test(dataProvider = "CompareDifferentDuplicatesProvider", expectedExceptions = UserException.class)
    public void differentBamTest(File firstBam, File secondBam) throws Exception {
        // These files are the same and should produce no diffs.
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(firstBam.getCanonicalPath());
        args.add("--" + "I2");
        args.add(secondBam.getCanonicalPath());
        args.add("--" + "throwOnDiff");
        args.add("true");

        this.runCommandLine(args.getArgsArray());
    }
}