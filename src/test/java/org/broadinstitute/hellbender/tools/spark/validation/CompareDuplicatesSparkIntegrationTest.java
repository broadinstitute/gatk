package org.broadinstitute.hellbender.tools.spark.validation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.ValidationStringency;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadQueryNameComparator;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.ReadTestUtils;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

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

    @Test(dataProvider = "CompareIdenticalDuplicatesProvider", groups = "spark")
    public void identicalBamTest(File firstBam, File secondBam) throws Exception {
        // These files are the same and should produce no diffs.

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(firstBam.getCanonicalPath());
        args.add("--" + CompareDuplicatesSpark.INPUT_2_SHORT_NAME);
        args.add(secondBam.getCanonicalPath());
        args.add("--" + CompareDuplicatesSpark.THROW_ON_DIFF_LONG_NAME);
        args.add("true");

        this.runCommandLine(args.getArgsArray());
    }

    @Test(dataProvider = "CompareDifferentDuplicatesProvider", expectedExceptions = UserException.class, groups = "spark")
    public void differentBamTest(File firstBam, File secondBam) throws Exception {
        // These files are not the same and should throw an exception
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(firstBam.getCanonicalPath());
        args.add("--" + CompareDuplicatesSpark.INPUT_2_SHORT_NAME);
        args.add(secondBam.getCanonicalPath());
        args.add("--" + CompareDuplicatesSpark.THROW_ON_DIFF_LONG_NAME);
        args.add("true");

        this.runCommandLine(args.getArgsArray());
    }

    @Test( groups = "spark")
    public void testOutputFile() throws Exception {
        final String resourceDir = getTestDataDir() + "/validation/";
        final File gatk4Bam = new File(resourceDir, "tmp.gatk4.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");
        final File dupfreeBam = new File(resourceDir, "unDuplicateMarked.gatk4.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam");
        final File validation1 = new File(resourceDir, "compareDuplicatesOutputTest.output1.bam");
        final File validation2 = new File(resourceDir, "compareDuplicatesOutputTest.output2.bam");
        final File diffbam1 = createTempFile("differentBamTest1", ".bam");
        final File diffbam2 = createTempFile("differentBamTest2", ".bam");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(gatk4Bam.getCanonicalPath());
        args.add("--" + CompareDuplicatesSpark.INPUT_2_SHORT_NAME);
        args.add(dupfreeBam.getCanonicalPath());
        args.add("-O");
        args.add(diffbam1);
        args.add("-O2");
        args.add(diffbam2);

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(diffbam1, validation1,  ValidationStringency.SILENT );
        SamAssertionUtils.assertSamsEqual(diffbam2, validation2,  ValidationStringency.SILENT );

        Pair<SAMFileHeader, List<GATKRead>> bam1Diffs = ReadTestUtils.readEntireBamIntoMemory(diffbam1.getPath());
        Pair<SAMFileHeader, List<GATKRead>> bam2Diffs = ReadTestUtils.readEntireBamIntoMemory(diffbam2.getPath());

        bam1Diffs.getRight().sort(new ReadQueryNameComparator());
        bam2Diffs.getRight().sort(new ReadQueryNameComparator());

        Assert.assertEquals(bam1Diffs.getRight().size(), bam2Diffs.getRight().size());
        for (int i = 0; i < bam1Diffs.getRight().size(); i++) {
            GATKRead left = bam1Diffs.getRight().get(i);
            GATKRead right = bam2Diffs.getRight().get(i);
            Assert.assertEquals(left.getName(), right.getName());
            Assert.assertEquals(left.getStart(), right.getStart());
            Assert.assertFalse(right.isDuplicate());
        }
    }
}