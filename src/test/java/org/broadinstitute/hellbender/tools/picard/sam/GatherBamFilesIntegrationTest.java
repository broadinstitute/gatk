package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class GatherBamFilesIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/GatherBamFiles");

    @Override
    public String getTestedClassName() {
        return GatherBamFiles.class.getSimpleName();
    }

    @DataProvider(name="gathering")
    public Object[][] filterTestData() {
        final File origBam = new File(TEST_DATA_DIR, "orig.bam");
        final List<File> splitBamsUnmappedLast = Arrays.asList(
                new File(TEST_DATA_DIR, "indchr1.bam"),
                new File(TEST_DATA_DIR, "indchr2.bam"),
                new File(TEST_DATA_DIR, "indchr3.bam"),
                new File(TEST_DATA_DIR, "indchr4.bam"),
                new File(TEST_DATA_DIR, "indchr5.bam"),
                new File(TEST_DATA_DIR, "indchr6.bam"),
                new File(TEST_DATA_DIR, "indchr7.bam"),
                new File(TEST_DATA_DIR, "indchr8.bam"),
                new File(TEST_DATA_DIR, "indUnknownChrom.bam")
        );

        final List<File> splitBamsUnmappedFirst = Arrays.asList(
                new File(TEST_DATA_DIR, "indUnknownChrom.bam"),
                new File(TEST_DATA_DIR, "indchr1.bam"),
                new File(TEST_DATA_DIR, "indchr2.bam"),
                new File(TEST_DATA_DIR, "indchr3.bam"),
                new File(TEST_DATA_DIR, "indchr4.bam"),
                new File(TEST_DATA_DIR, "indchr5.bam"),
                new File(TEST_DATA_DIR, "indchr6.bam"),
                new File(TEST_DATA_DIR, "indchr7.bam"),
                new File(TEST_DATA_DIR, "indchr8.bam")
        );
        final List<File> splitBamsUnmappedMissing = Arrays.asList(
                new File(TEST_DATA_DIR, "indchr1.bam"),
                new File(TEST_DATA_DIR, "indchr2.bam"),
                new File(TEST_DATA_DIR, "indchr3.bam"),
                new File(TEST_DATA_DIR, "indchr4.bam"),
                new File(TEST_DATA_DIR, "indchr5.bam"),
                new File(TEST_DATA_DIR, "indchr6.bam"),
                new File(TEST_DATA_DIR, "indchr7.bam"),
                new File(TEST_DATA_DIR, "indchr8.bam")
        );

        return new Object[][]{
                {splitBamsUnmappedLast,  origBam,   true, true},
                {splitBamsUnmappedFirst, origBam,   false, true},
                {splitBamsUnmappedMissing, origBam, false, false},
        };
    }

    @Test(dataProvider = "gathering")
    public void testTheGathering(final List<File> bams, final File original, final boolean expectEqual, final boolean expectEqualLenient) throws IOException {
        final File outputFile = BaseTest.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        final List<String> args = new ArrayList<>();
        for (final File splitBam : bams) {
            args.add("--INPUT");
            args.add(splitBam.getAbsolutePath());
        }
        args.add("--OUTPUT");
        args.add(outputFile.getAbsolutePath());
        runCommandLine(args);
        if (expectEqual) {
            Assert.assertNull(SamAssertionUtils.samsEqualStringent(original, outputFile, ValidationStringency.DEFAULT_STRINGENCY, null));
        } else {
            Assert.assertNotNull(SamAssertionUtils.samsEqualStringent(original, outputFile, ValidationStringency.DEFAULT_STRINGENCY, null));
        }
        if (expectEqualLenient) {
            Assert.assertNull(SamAssertionUtils.samsEqualLenient(original, outputFile, ValidationStringency.DEFAULT_STRINGENCY, null));
        } else {
            Assert.assertNotNull(SamAssertionUtils.samsEqualLenient(original, outputFile, ValidationStringency.DEFAULT_STRINGENCY, null));
        }
    }
}
