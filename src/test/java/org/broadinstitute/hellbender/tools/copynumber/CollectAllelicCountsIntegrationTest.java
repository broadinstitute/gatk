package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

/**
 * Integration test for {@link CollectAllelicCounts}.  Uses a BAM with sites generated from hg19mini using wgsim.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CollectAllelicCountsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_SUB_DIR = toolsTestDir + "copynumber/allelic";
    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-normal.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-tumor.bam");
    private static final File SITES_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-sites.interval_list");
    private static final File REFERENCE_FILE = new File(hg19MiniReference);
    private static final String NORMAL_SAMPLE_NAME_EXPECTED = "20";
    private static final String TUMOR_SAMPLE_NAME_EXPECTED = "20";

    @DataProvider(name = "testData")
    public Object[][] testData() {
        //counts from IGV with minMQ = 30 and minBQ = 20
        final AllelicCountCollection normalCountsExpected = new AllelicCountCollection(
                new SimpleSampleMetadata(NORMAL_SAMPLE_NAME_EXPECTED),
                Arrays.asList(
                        new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0, Nucleotide.G, Nucleotide.N),
                        new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4, Nucleotide.G, Nucleotide.A),
                        new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6, Nucleotide.G, Nucleotide.T),
                        new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18, Nucleotide.T, Nucleotide.C),
                        new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8, Nucleotide.C, Nucleotide.A),
                        new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8, Nucleotide.T, Nucleotide.G),
                        new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4, Nucleotide.C, Nucleotide.A),
                        new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9, Nucleotide.T, Nucleotide.G),
                        new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5, Nucleotide.G, Nucleotide.C),
                        new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0, Nucleotide.G, Nucleotide.N),
                        new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3, Nucleotide.T, Nucleotide.A)));

        final AllelicCountCollection tumorCountsExpected = new AllelicCountCollection(
                new SimpleSampleMetadata(TUMOR_SAMPLE_NAME_EXPECTED),
                Arrays.asList(
                        new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0, Nucleotide.G, Nucleotide.N),
                        new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4, Nucleotide.G, Nucleotide.A),
                        new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6, Nucleotide.G, Nucleotide.T),
                        new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 17, Nucleotide.T, Nucleotide.C),
                        new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8, Nucleotide.C, Nucleotide.A),
                        new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8, Nucleotide.T, Nucleotide.G),
                        new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 3, Nucleotide.C, Nucleotide.A),
                        new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9, Nucleotide.T, Nucleotide.G),
                        new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5, Nucleotide.G, Nucleotide.C),
                        new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0, Nucleotide.G, Nucleotide.N),
                        new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3, Nucleotide.T, Nucleotide.A)));

        //counts from IGV with minMQ = 30 and minBQ = 20, without nucleotides
        final AllelicCountCollection normalCountsExpectedWithoutNucleotides = new AllelicCountCollection(
                new SimpleSampleMetadata(NORMAL_SAMPLE_NAME_EXPECTED),
                Arrays.asList(
                        new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0),
                        new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4),
                        new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6),
                        new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18),
                        new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8),
                        new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8),
                        new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4),
                        new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9),
                        new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5),
                        new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0),
                        new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3)));

        final AllelicCountCollection tumorCountsExpectedWithoutNucleotides = new AllelicCountCollection(
                new SimpleSampleMetadata(TUMOR_SAMPLE_NAME_EXPECTED),
                Arrays.asList(
                        new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0),
                        new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4),
                        new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6),
                        new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 17),
                        new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8),
                        new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8),
                        new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 3),
                        new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9),
                        new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5),
                        new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0),
                        new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3)));

        return new Object[][]{
                {NORMAL_BAM_FILE, normalCountsExpected},
                {TUMOR_BAM_FILE, tumorCountsExpected},
                {NORMAL_BAM_FILE, normalCountsExpectedWithoutNucleotides},
                {TUMOR_BAM_FILE, tumorCountsExpectedWithoutNucleotides}
        };
    }

    @Test(dataProvider = "testData")
    public void test(final File inputBAMFile,
                     final AllelicCountCollection countsExpected) {
        final File outputFile = createTempFile("collect-allelic-counts-test-output", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputBAMFile.getAbsolutePath(),
                "-L", SITES_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        final AllelicCountCollection countsResult = new AllelicCountCollection(outputFile);
        Assert.assertEquals(countsExpected, countsResult);
        Assert.assertNotSame(countsExpected, countsResult);
    }
}