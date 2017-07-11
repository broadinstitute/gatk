package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Integration test for {@link CollectAllelicCounts}.  Uses BAM and SNP files generated from hg19mini using wgsim.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CollectAllelicCountsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/allelic";
    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-normal.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-tumor.bam");
    private static final File NON_STRICT_BAM_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-simple-overhang.sam");
    private static final File SITES_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-sites.interval_list");
    private static final File REF_FILE = new File(hg19MiniReference);

    @DataProvider(name = "testData")
    public Object[][] testData() throws IOException {
        //counts from IGV with minMQ = 30 and minBQ = 20
        final AllelicCountCollection normalCountsExpected = new AllelicCountCollection();
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3));

        final AllelicCountCollection tumorCountsExpected = new AllelicCountCollection();
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 17));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 3));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3));

        return new Object[][]{
                {NORMAL_BAM_FILE, normalCountsExpected},
                {TUMOR_BAM_FILE, tumorCountsExpected}
        };
    }

    @Test(dataProvider = "testData")
    public void test(final File inputBAMFile,
                     final AllelicCountCollection countsExpected) {
        final File outputFile = createTempFile("collect-allelic-counts-test-output", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputBAMFile.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SITES_FILE_SHORT_NAME, SITES_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        final AllelicCountCollection countsResult = new AllelicCountCollection(outputFile);
        Assert.assertEquals(countsExpected, countsResult);
    }

    //Regression test for https://github.com/broadinstitute/gatk-protected/issues/373
    @Test(expectedExceptions = UserException.class)
    public void testNonStrictBAM() {
        final File outputFile = createTempFile("collect-allelic-counts-test-output", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, NON_STRICT_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SITES_FILE_SHORT_NAME, SITES_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_SHORT_NAME, ValidationStringency.STRICT.toString()
        };
        runCommandLine(arguments);
        //should catch SAMFormatException and throw new UserException with --readValidationStringency STRICT
    }

    //Regression test for https://github.com/broadinstitute/gatk-protected/issues/373
    @Test
    public void testNonStrictBAMWithSilentValidationStringency() {
        final File outputFile = createTempFile("collect-allelic-counts-test-output", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, NON_STRICT_BAM_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SITES_FILE_SHORT_NAME, SITES_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_SHORT_NAME, ValidationStringency.SILENT.toString()
        };
        runCommandLine(arguments);
        //should complete successfully with --readValidationStringency SILENT
    }
}