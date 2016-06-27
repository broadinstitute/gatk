package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.pulldown.Pulldown;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Integration test for {@link GetBayesianHetCoverage}. Uses BAM and SNP files generated from hg19mini using wgsim.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GetBayesianHetCoverageIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";

    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "normal.sorted.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR, "tumor.sorted.bam");
    private static final File SNP_FILE = new File(TEST_SUB_DIR, "common_SNP.interval_list");
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

    @Test(expectedExceptions = UserException.class)
    public void testBadJobArgs_1() {
        final File normalOutputFile = createTempFile("normal-test", ".tsv");

        /* bad job: no normal BAM file */
        final String[] argumentsNoBam = {
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME, normalOutputFile.getAbsolutePath()
        };
        runCommandLine(argumentsNoBam);
    }

    @Test(expectedExceptions = UserException.class)
    public void testBadJobArgs_2() {
        /* bad job: no normal het pulldown output file */
        final String[] argumentsNormalNoHetOutput = {
                "-" + ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_SHORT_NAME, NORMAL_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath()
        };
        runCommandLine(argumentsNormalNoHetOutput);
    }

    @Test(expectedExceptions = UserException.class)
    public void testBadJobArgs_3() {
        /* bad job: no tumor het pulldown output file */
        final String[] argumentsTumorNoHetOutput = {
                "-" + ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME, NORMAL_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath()
        };
        runCommandLine(argumentsTumorNoHetOutput);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadPriorParameters_1() {
        final File tumorOutputFile = createTempFile("tumor-test", ".tsv");

        /* bad choice of prior parameters: max abnormal fraction > 1 */
        final String[] argumentsBadPrior = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.MAXIMUM_ABNORMAL_FRACTION_SHORT_NAME, Double.toString(1.1)
        };
        runCommandLine(argumentsBadPrior);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadPriorParameters_2() {
        final File tumorOutputFile = createTempFile("tumor-test", ".tsv");

        /* bad choice of prior parameters: min abnormal fraction < 0 */
        final String[] argumentsBadPrior = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.MINIMUM_ABNORMAL_FRACTION_SHORT_NAME, Double.toString(-0.1)
        };
        runCommandLine(argumentsBadPrior);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadPriorParameters_3() {
        final File tumorOutputFile = createTempFile("tumor-test", ".tsv");

        /* bad choice of prior parameters: max abnormal fraction < min abnormal fraction */
        final String[] argumentsBadPrior = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.MINIMUM_ABNORMAL_FRACTION_SHORT_NAME, Double.toString(0.5),
                "-" + GetBayesianHetCoverage.MAXIMUM_ABNORMAL_FRACTION_SHORT_NAME, Double.toString(0.4)
        };
        runCommandLine(argumentsBadPrior);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadPriorParameters_4() {
        final File tumorOutputFile = createTempFile("tumor-test", ".tsv");

        /* bad choice of prior parameters: maximum copy number < 0 */
        final String[] argumentsBadPrior = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.MAXIMUM_COPY_NUMBER_SHORT_NAME, Integer.toString(-1)
        };
        runCommandLine(argumentsBadPrior);
    }

    @Test
    public void testNormalJob() {
        final File normalOutputFile = createTempFile("normal-test", ".tsv");
        Pulldown normalHetPulldownExpected, normalHetPulldownResult;

        /* test 1: tight calling stringency */
        final String[] args_1 = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_SHORT_NAME, NORMAL_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME, normalOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.READ_DEPTH_THRESHOLD_SHORT_NAME, Integer.toString(10),
                "-" + GetBayesianHetCoverage.HET_CALLING_STRINGENCY_SHORT_NAME, Double.toString(15.0)
        };
        runCommandLine(args_1);

        normalHetPulldownExpected = new Pulldown(normalHeader);
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                Nucleotide.T, Nucleotide.G, 17, 39.99));

        normalHetPulldownResult = new Pulldown(normalOutputFile, normalHeader);

        Assert.assertEquals(normalHetPulldownExpected, normalHetPulldownResult);

        /* test 2: loose calling stringency */
        final String[] args_2 = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_SHORT_NAME, NORMAL_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME, normalOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.READ_DEPTH_THRESHOLD_SHORT_NAME, Integer.toString(10),
                "-" + GetBayesianHetCoverage.HET_CALLING_STRINGENCY_SHORT_NAME, Double.toString(2.0)
        };
        runCommandLine(args_2);

        normalHetPulldownExpected = new Pulldown(normalHeader);
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4,
                Nucleotide.G, Nucleotide.A, 11, 18.60));
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6,
                Nucleotide.G, Nucleotide.T, 14, 29.29));
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                Nucleotide.T, Nucleotide.G, 17, 39.99));
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9,
                Nucleotide.T, Nucleotide.G, 15, 28.60));
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5,
                Nucleotide.G, Nucleotide.C, 11, 24.99));
        normalHetPulldownResult = new Pulldown(normalOutputFile, normalHeader);

        Assert.assertEquals(normalHetPulldownExpected, normalHetPulldownResult);
    }

    @Test
    public void testTumorJob() {
        final File tumorOutputFile = createTempFile("tumor-test", ".tsv");
        Pulldown tumorHetPulldownExpected, tumorHetPulldownResult;

        /* test 1: tight calling stringency */
        final String[] args_1 = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.READ_DEPTH_THRESHOLD_SHORT_NAME, Integer.toString(10),
                "-" + GetBayesianHetCoverage.HET_CALLING_STRINGENCY_SHORT_NAME, Double.toString(10.0),
        };
        runCommandLine(args_1);

        tumorHetPulldownExpected = new Pulldown(tumorHeader);
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6,
                Nucleotide.G, Nucleotide.T, 14, 24.85));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                Nucleotide.T, Nucleotide.G, 17, 34.03));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9,
                Nucleotide.T, Nucleotide.G, 15, 24.23));

        tumorHetPulldownResult = new Pulldown(tumorOutputFile, tumorHeader);

        Assert.assertEquals(tumorHetPulldownExpected, tumorHetPulldownResult);

        /* test 2: tight calling stringency */
        final String[] args_2 = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.READ_DEPTH_THRESHOLD_SHORT_NAME, Integer.toString(10),
                "-" + GetBayesianHetCoverage.HET_CALLING_STRINGENCY_SHORT_NAME, Double.toString(5.0)
        };
        runCommandLine(args_2);

        tumorHetPulldownExpected = new Pulldown(tumorHeader);
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4,
                Nucleotide.G, Nucleotide.A, 11, 15.72));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6,
                Nucleotide.G, Nucleotide.T, 14, 24.85));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                Nucleotide.T, Nucleotide.G, 17, 34.03));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9,
                Nucleotide.T, Nucleotide.G, 15, 24.23));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5,
                Nucleotide.G, Nucleotide.C, 11, 21.23));

        tumorHetPulldownResult = new Pulldown(tumorOutputFile, tumorHeader);

        Assert.assertEquals(tumorHetPulldownExpected, tumorHetPulldownResult);
    }

    @Test
    public void testMatchedNormalTumorJob() {
        final File normalOutputFile = createTempFile("normal-test", ".tsv");
        final File tumorOutputFile = createTempFile("tumor-test", ".tsv");
        Pulldown tumorHetPulldownExpected, tumorHetPulldownResult;
        Pulldown normalHetPulldownExpected, normalHetPulldownResult;

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME, SNP_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.NORMAL_BAM_FILE_SHORT_NAME, NORMAL_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME, normalOutputFile.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_BAM_FILE_SHORT_NAME, TUMOR_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, tumorOutputFile.getAbsolutePath(),
                "-" + GetBayesianHetCoverage.READ_DEPTH_THRESHOLD_SHORT_NAME, Integer.toString(10),
                "-" + GetBayesianHetCoverage.HET_CALLING_STRINGENCY_SHORT_NAME, Double.toString(10.0)
        };
        runCommandLine(arguments);

        normalHetPulldownResult = new Pulldown(normalOutputFile, normalHeader);
        tumorHetPulldownResult = new Pulldown(tumorOutputFile, tumorHeader);

        normalHetPulldownExpected = new Pulldown(normalHeader);
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6,
                Nucleotide.G, Nucleotide.T, 14, 29.29));
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                Nucleotide.T, Nucleotide.G, 17, 39.98));
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9,
                Nucleotide.T, Nucleotide.G, 15, 28.60));
        normalHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5,
                Nucleotide.G, Nucleotide.C, 11, 24.99));

        tumorHetPulldownExpected = new Pulldown(tumorHeader);
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6,
                Nucleotide.G, Nucleotide.T, 14));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                Nucleotide.T, Nucleotide.G, 17));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9,
                Nucleotide.T, Nucleotide.G, 15));
        tumorHetPulldownExpected.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5,
                Nucleotide.G, Nucleotide.C, 11));

        Assert.assertEquals(normalHetPulldownExpected, normalHetPulldownResult);
        Assert.assertEquals(tumorHetPulldownExpected, tumorHetPulldownResult);
    }
}
