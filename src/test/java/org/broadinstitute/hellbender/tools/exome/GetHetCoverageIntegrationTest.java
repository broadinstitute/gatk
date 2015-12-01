package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException.CommandLineException;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Integration test for {@link GetHetCoverage}.  Uses BAM and SNP files generated from hg19mini using wgsim.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GetHetCoverageIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR + "normal.sorted.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR + "tumor.sorted.bam");
    private static final File SNP_FILE = new File(TEST_SUB_DIR + "common_SNP.interval_list");
    private static final File REF_FILE = new File(hg19MiniReference);

    private static SAMFileHeader normalHeader;
    private static SAMFileHeader tumorHeader;

    @BeforeClass
    public void initHeaders() throws IOException {
        try (final SamReader normalBamReader = SamReaderFactory.makeDefault().open(NORMAL_BAM_FILE);
             final SamReader tumorBamReader = SamReaderFactory.makeDefault().open(TUMOR_BAM_FILE)) {
            normalHeader = normalBamReader.getFileHeader();
            tumorHeader = tumorBamReader.getFileHeader();
        }
    }

    @Test
    public void testGetHetCoverage() {
        final File normalOutputFile = createTempFile("normal-test", ".txt");
        final File tumorOutputFile = createTempFile("tumor-test", ".txt");

        final String[] arguments = {
                "-" + GetHetCoverage.NORMAL_BAM_FILE_SHORT_NAME, NORMAL_BAM_FILE.getAbsolutePath(),
                "-" + GetHetCoverage.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + GetHetCoverage.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + GetHetCoverage.NORMAL_HET_REF_ALT_COUNTS_FILE_SHORT_NAME, normalOutputFile.getAbsolutePath(),
                "-" + GetHetCoverage.TUMOR_HET_REF_ALT_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        final Pulldown normalOutputPulldownCLP = new Pulldown(normalOutputFile, normalHeader);
        final Pulldown tumorOutputPulldownCLP = new Pulldown(tumorOutputFile, tumorHeader);

        Pulldown normalHetPulldown = new Pulldown(normalHeader);
        normalHetPulldown.add(new Interval("1", 10736, 10736), 9, 2);
        normalHetPulldown.add(new Interval("1", 11522, 11522), 7, 4);
        normalHetPulldown.add(new Interval("1", 12098, 12098), 8, 6);
        normalHetPulldown.add(new Interval("1", 14630, 14630), 9, 8);
        normalHetPulldown.add(new Interval("2", 14689, 14689), 6, 9);
        normalHetPulldown.add(new Interval("2", 14982, 14982), 6, 5);

        Pulldown tumorHetPulldown = new Pulldown(tumorHeader);
        tumorHetPulldown.add(new Interval("1", 11522, 11522), 7, 4);
        tumorHetPulldown.add(new Interval("1", 12098, 12098), 8, 6);
        tumorHetPulldown.add(new Interval("1", 14630, 14630), 9, 8);
        tumorHetPulldown.add(new Interval("2", 14689, 14689), 6, 9);
        tumorHetPulldown.add(new Interval("2", 14982, 14982), 6, 5);

        Assert.assertEquals(normalHetPulldown, normalOutputPulldownCLP);
        Assert.assertEquals(tumorHetPulldown, tumorOutputPulldownCLP);
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testBadpValue() {
        final File normalOutputFile = createTempFile("normal-test",".txt");
        final File tumorOutputFile = createTempFile("tumor-test",".txt");

        final String[] arguments = {
                "-" + GetHetCoverage.NORMAL_BAM_FILE_SHORT_NAME, NORMAL_BAM_FILE.getAbsolutePath(),
                "-" + GetHetCoverage.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + GetHetCoverage.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + GetHetCoverage.NORMAL_HET_REF_ALT_COUNTS_FILE_SHORT_NAME, normalOutputFile.getAbsolutePath(),
                "-" + GetHetCoverage.TUMOR_HET_REF_ALT_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath(),
                "-" + GetHetCoverage.PVALUE_THRESHOLD_SHORT_NAME, Double.toString(-10)
        };
        runCommandLine(arguments);
    }
}
